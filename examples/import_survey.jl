using MagneticGeophysics

length(ARGS) >= 4 || error("usage: julia --project=. examples/import_survey.jl survey.csv F I D [mt_data.dat [mt_model.rho]]")
survey_file  = ARGS[1]
field        = Tuple(parse.(Float64, ARGS[2:4]))      # inducing field F (nT), I, D (degrees, D from true north)
mt_data_file = length(ARGS) >= 5 ? ARGS[5] : ""
mt_model     = length(ARGS) >= 6 ? ARGS[6] : ""
out_file     = splitext(survey_file)[1] * ".dat"

component      = :tmi              # :tmi, :bx, :by, :bz, :amplitude, :rtp, :dtmi_dz, :dbz_dz
crs            = "EPSG:4326"       # of easting/northing columns
sensor_azimuth = 0.0               # instrument x axis of Bx/By data, clockwise from north
# edit to the column names of the survey file (lat/lon, easting/northing or north/east)
columns = (lat = "LAT", lon = "LON", elevation = "ALT", value = "TMI", error = "ERR", line = "LINE")

frame = isempty(mt_data_file) ? nothing : load_frame(mt_data_file)
d = load_data_mag(survey_file; component, field, columns, crs, frame, sensor_azimuth)
if !isempty(mt_model)
    m = load_mesh(mt_model)
    m.rotation == 0 || (d = rotate_data(d, m.rotation))
    check_frame(d, m)
end
println(write_data_mag(out_file, d; comment = "imported from $(basename(survey_file))"))
