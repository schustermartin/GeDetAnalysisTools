# GeDetAnalysisTools.jl 

One has to create two environement variables:

- `GEDET_USER_DATA_PATH`: Path to the data files. (e.g.: `export GEDET_USER_DATA_PATH="/remote/ceph/group/gedet/data/lab"`) 
Inside this path their must be the correct data structure (NOT_YET_DOCUMENTET).
- `GEDET_USER_OUTPUT_PATH` Path to where output should be produced. (e.g.: `export GEDET_USER_OUTPUT_PATH="/remote/ceph/user/<u>/<username>/analysis_directory"`)

## Example usage of data file handling:
```
julia> using Plots; pyplot(); using GeDetAnalysisTools

julia> data_set_name = "2018-10-11_ea2cd342_lm_GALATEA_SuSie_rotational_scan_tsrc_Am241_ssrc_Am241_r23_z35";

julia> dataset = data_set_from_conv_data(data_set_name)

2018-10-11_ea2cd342_lm_GALATEA_SuSie_rotational_scan_tsrc_Am241_ssrc_Am241_r23_z35
contains 37 measurements

julia> m = dataset[1]

name: R_023.001mm_Z_034.998mm_Phi_004.999deg_T_318.57K_p_8.6e-06mbar_t_3600sec
data_set_name: 2018-10-11_ea2cd342_lm_GALATEA_SuSie_rotational_scan_tsrc_Am241_ssrc_Am241_r23_z35
path_to_raw_data: /remote/ceph/group/gedet/data/lab/2018/2018-10-11_ea2cd342_lm_GALATEA_SuSie_rotational_scan_tsrc_Am241_ssrc_Am241_r23_z35/raw_data
results_path: /remote/ceph/user/l/lhauert/analysis_directory/2018-10-11_ea2cd342_lm_GALATEA_SuSie_rotational_scan_tsrc_Am241_ssrc_Am241_r23_z35
datetime: 2018-10-11T14:56:19
side_source: Am241
top_source: Am241
outside_source: unknown
motor_pos_z: 34.998
motor_pos_r: 23.001
motor_pos_phi: 4.999
sidesourcefacingsegment: -1
topsourcefacingsegment: -1
detector: Detector("noname", 2, [1, 2], [1, 2], Int64[], Int64[], Int64[], Int64[])
measurement_duration: 3600
daq: DAQ("noname", -1, -1, NaN, -1, -1)
daq_livetime: -1.0
temperature: 318.57
pressure: 8.6e-6
tau_decay_constants: Float64[]
new_data_structure: true

julia> daq_cal_energies = GeDetAnalysisTools.get_quick_calibrated_daq_energies(m);

julia> plot_energy_histograms(daq_cal_energies, SUSIE(), edges=0:1:6000) # @userplot recipe
```
