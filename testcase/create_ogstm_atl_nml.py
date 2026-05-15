import os
import pyfabm

CODEPATH = '../../'
CODEPATH = CODEPATH.replace("~", os.getenv("HOME"))
fabm_yaml = "/leonardo_work/OGS_test2528_0/plazzari/OGSTM-FABM/ModelBuild/ogstm/testcase/TEST02/wrkdir/MODEL/fabm.yaml"


def generate_atl_nml(fabm_yaml_path, output_file="atl.nml"):
    model = pyfabm.Model(fabm_yaml_path)

    state_vars = [v.name.replace('/', '_') for v in model.state_variables]
    n_vars = len(state_vars)

    with open(output_file, 'w') as f:
        f.write("&VARS_DIMENSION\n")
        f.write(f"    n_vars = {n_vars}\n")
        f.write("/\n")
        f.write("\n")
        f.write("&CORE\n")
        f.write("\n")
        for i, name in enumerate(state_vars, start=1):
            f.write(f'    vars({i}) = "{name}"\n')
        f.write("/\n")

    print(f"Generated {output_file} with {n_vars} state variables.")


if __name__ == "__main__":
    generate_atl_nml(fabm_yaml)
