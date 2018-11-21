"""

"""

import swmm_output as so


link_id = "C7"
 
# Initialize Output Object Handle
sm_handle = so.smo_init()

# Binary File
out_file = r'./Set1.out'

# Open Binary File
so.smo_open(sm_handle, out_file)

# Get Periods
nPeriods = so.smo_get_times(sm_handle, so.Time.NUM_PERIODS)
# Get Timestep (seconds)
tstep = so.smo_get_times(sm_handle, so.Time.REPORT_STEP)

# Search for Index of the object in the binary file.
link_dict = {}
num_links = so.smo_get_project_size(sm_handle)[so.ElementType.LINK.value]

for i in range(num_links):
    lnkname = so.smo_get_element_name(sm_handle, so.ElementType.LINK, i)
    link_dict[lnkname] = i


for i, val in enumerate(so.smo_get_link_series(sm_handle, link_dict[link_id],
                                               so.LinkAttribute.FLOW_RATE,
                                               0, nPeriods)):
    print(val)

# close output file
sm_handle = so.smo_close()



