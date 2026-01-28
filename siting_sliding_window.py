import os
import arcpy
from arcpy.sa import (
    Raster, FocalStatistics, NbrCircle, SetNull, Con, IsNull, CreateConstantRaster
)

# =============================================================================
# Sliding-window siting with iterative biomass-competition deduction (ArcGIS Pro 3.x)
# Outputs:
#   - outputs/focal/focal_sum_0001.tif ... focal_sum_NNNN.tif
#   - outputs/focal/topXXXX_sites.gdb/Sites_Top*_Points (point FC)
#   - outputs/focal/topXXXX_sites.gdb/Sites_Top*_Buffers (polygon FC)
# Notes:
#   - Designed to run in the ArcGIS Pro Python window.
#   - Keeps intermediate rasters in a FileGDB to reduce TIFF locking issues.
# =============================================================================

# -------------------------
# User parameters (edit here)
# -------------------------
IN_RASTER_REL = os.path.join("data", "biomethane_map.tif")
OUT_DIR_REL   = os.path.join("outputs", "focal")
N_SITES       = 1500
RADIUS_CELLS  = 25
# Set to True only for debugging (e.g., to inspect hotspots or troubleshoot ArcGIS errors); note this will generate many large intermediate focal_sum_####.tif files—disable after debugging and delete intermediates if needed.
SAVE_FOCAL_TIFS = False  
# -------------------------

def _project_root():
    # In ArcGIS Pro Python window, __file__ is not defined; fall back to CWD.
    return os.path.abspath(os.path.dirname(__file__)) if "__file__" in globals() else os.getcwd()

def _safe_delete(path):
    if arcpy.Exists(path):
        try:
            arcpy.management.Delete(path)
        except Exception:
            pass

def run():
    project_root = _project_root()
    in_raster = os.path.join(project_root, IN_RASTER_REL)
    out_dir   = os.path.join(project_root, OUT_DIR_REL)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    arcpy.CheckOutExtension("Spatial")
    arcpy.env.overwriteOutput = True
    arcpy.env.addOutputsToMap = False

    # Avoid occasional disk/CopyRaster issues by disabling parallel processing
    arcpy.env.parallelProcessingFactor = "0%"

    # Align all derived rasters with the input raster
    arcpy.env.snapRaster = in_raster
    arcpy.env.cellSize   = in_raster
    arcpy.env.extent     = in_raster

    # Store working rasters and outputs in a FileGDB (more stable than many TIFF writes)
    out_gdb_name = f"top{N_SITES}_sites.gdb"
    out_gdb = os.path.join(out_dir, out_gdb_name)
    _safe_delete(out_gdb)
    arcpy.management.CreateFileGDB(out_dir, out_gdb_name)
    arcpy.env.workspace = out_gdb
    arcpy.env.scratchWorkspace = out_gdb

    base0 = Raster(in_raster)
    cell_size = float(base0.meanCellWidth)
    radius_deg = RADIUS_CELLS * cell_size
    sr = arcpy.Describe(in_raster).spatialReference

    print(f"Input raster: {in_raster}")
    print(f"CellSize = {cell_size} | RadiusCells = {RADIUS_CELLS} | RadiusDeg = {radius_deg}")
    print(f"SpatialReference type: {sr.type}")

    # Feature classes for selected sites
    sites_fc   = os.path.join(out_gdb, f"Sites_Top{N_SITES}_Points")
    buffers_fc = os.path.join(out_gdb, f"Sites_Top{N_SITES}_Buffers")

    arcpy.management.CreateFeatureclass(out_gdb, os.path.basename(sites_fc), "POINT", spatial_reference=sr)
    arcpy.management.AddField(sites_fc, "Iter", "LONG")
    arcpy.management.AddField(sites_fc, "MaxSum", "DOUBLE")

    arcpy.management.CreateFeatureclass(out_gdb, os.path.basename(buffers_fc), "POLYGON", spatial_reference=sr)
    arcpy.management.AddField(buffers_fc, "Iter", "LONG")
    arcpy.management.AddField(buffers_fc, "MaxSum", "DOUBLE")

    # Working rasters (toggle A/B)
    work_a = os.path.join(out_gdb, "working_a")
    work_b = os.path.join(out_gdb, "working_b")
    excl_a = os.path.join(out_gdb, "exclude_a")
    excl_b = os.path.join(out_gdb, "exclude_b")

    # Initialize working raster
    try:
        arcpy.management.Copy(in_raster, work_a)
    except Exception:
        Raster(in_raster).save(work_a)

    # Initialize exclusion mask (0 everywhere)
    CreateConstantRaster(0, "INTEGER", cell_size).save(excl_a)

    neighborhood = NbrCircle(RADIUS_CELLS, "CELL")

    # Temporary datasets (overwritten each iteration)
    tmp_search   = os.path.join(out_gdb, "tmp_focal_search")
    tmp_maxcells = os.path.join(out_gdb, "tmp_maxcells")
    tmp_mask     = os.path.join(out_gdb, "tmp_mask")
    tmp_pts      = os.path.join(out_gdb, "tmp_max_pts")
    tmp_onept    = os.path.join(out_gdb, "tmp_onept")
    tmp_buf      = os.path.join(out_gdb, "tmp_buf")

    work_in, work_out = work_a, work_b
    excl_in, excl_out = excl_a, excl_b

    # Buffer distance string
    if sr.type == "Geographic":
        buffer_dist = f"{radius_deg} DecimalDegrees"
    else:
        buffer_dist = f"{RADIUS_CELLS * cell_size}"

    for k in range(1, N_SITES + 1):
        wr  = Raster(work_in)
        exr = Raster(excl_in)

        focal = FocalStatistics(wr, neighborhood, "SUM", "DATA")

        if SAVE_FOCAL_TIFS:
            focal_out = os.path.join(out_dir, f"focal_sum_{k:04d}.tif")
            focal.save(focal_out)

        # Exclude previously selected areas (exclusion mask == 1 -> NoData)
        focal_search = SetNull(exr, focal, "VALUE = 1")

        for p in [tmp_search, tmp_maxcells, tmp_mask, tmp_pts, tmp_onept, tmp_buf]:
            _safe_delete(p)

        focal_search.save(tmp_search)
        arcpy.management.CalculateStatistics(tmp_search)

        # Maximum value
        try:
            max_val = float(arcpy.management.GetRasterProperties(tmp_search, "MAXIMUM").getOutput(0))
        except Exception:
            print(f"[{k}] No valid MAXIMUM (all NoData). Stop.")
            break

        if max_val <= 0:
            print(f"[{k}] MAXIMUM <= 0 (max={max_val}). Stop.")
            break

        # Locate the (first) cell(s) reaching max_val -> point
        max_cells = Con(Raster(tmp_search) == max_val, 1)
        max_cells.save(tmp_maxcells)

        arcpy.conversion.RasterToPoint(tmp_maxcells, tmp_pts, "VALUE")

        xy = None
        with arcpy.da.SearchCursor(tmp_pts, ["SHAPE@XY"]) as cur:
            for row in cur:
                xy = row[0]
                break

        if xy is None:
            print(f"[{k}] Failed to find max location point. Stop.")
            break

        x, y = xy
        print(f"[{k}] MaxSum={max_val} at (x={x}, y={y})")

        pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)

        # Write point
        with arcpy.da.InsertCursor(sites_fc, ["SHAPE@", "Iter", "MaxSum"]) as icur:
            icur.insertRow([pt_geom, k, max_val])

        # Buffer
        arcpy.management.CopyFeatures([pt_geom], tmp_onept)
        arcpy.analysis.Buffer(tmp_onept, tmp_buf, buffer_dist, dissolve_option="NONE")

        # Write buffer polygon
        with arcpy.da.SearchCursor(tmp_buf, ["SHAPE@"]) as cur:
            for row in cur:
                with arcpy.da.InsertCursor(buffers_fc, ["SHAPE@", "Iter", "MaxSum"]) as icur:
                    icur.insertRow([row[0], k, max_val])
                break

        # Polygon -> raster -> binary mask
        arcpy.conversion.PolygonToRaster(
            in_features=tmp_buf,
            value_field="OBJECTID",
            out_rasterdataset=tmp_mask,
            cell_assignment="CELL_CENTER",
            cellsize=cell_size
        )
        mask_bin = Con(IsNull(Raster(tmp_mask)), 0, 1)

        # Deduct competition: set values inside buffer to 0
        updated_work = Con(mask_bin == 1, 0, wr)
        _safe_delete(work_out)
        updated_work.save(work_out)

        # Update exclusion mask: inside buffer = 1
        updated_excl = Con(mask_bin == 1, 1, exr)
        _safe_delete(excl_out)
        updated_excl.save(excl_out)

        # Swap for next iteration
        work_in, work_out = work_out, work_in
        excl_in, excl_out = excl_out, excl_in

        if k % 50 == 0:
            try:
                arcpy.management.ClearWorkspaceCache()
            except Exception:
                pass

    print("Done.")
    print("Points FC :", sites_fc)
    print("Buffers FC:", buffers_fc)
    if SAVE_FOCAL_TIFS:
        print(f"Focal TIFFs: {os.path.join(out_dir, 'focal_sum_0001.tif')} ...")

    arcpy.CheckInExtension("Spatial")


# Run
run()
