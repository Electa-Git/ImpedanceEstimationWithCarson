The yearly fifteen-minute-resolution end-user power profiles are taken from the Open Energy Data Initiative’s (OEDI) dataset [“End-Use Load Profiles for the U.S. Building Stock”](https://data.openei.org/submissions/4520).

The original OEDI's profiles only contain active power. To create reactive power, a random power factor between 0.95 and 0.998 is assigned for each user at each time step.

This repository only contains a subset of the original OEDI's dataset: representative profiles of 108 users from the State of New York are picked. Users are chosen that present average yearly active power between 0.058 kW and 2.8 kW, and maximum power not exceeding 12 kW. This creates a set in which some users have consumption similar to today's users without electric vehicles, heat pumps etc., whereas some users will have some of these devices. We note that impedance estimation quality tendentially increases with increasing power flows (as a consequence of increasing generation/demand). Thus, the purpose of this profile selection exercise is to "keep the exercise realistic" by filtering consumers with unrealistically/unusually high power (at least for today's standards), which would artifically improve the quality of impedance estimation.

The last part of the column names in our `profiles.csv` file, e.g., "_100106-2", matches the profile id from the OEDI database, so that it can be traced back to the original dataset.

Our `profiles.csv`, as well as the original OEDI dataset, have a CC BY 4.0 license.