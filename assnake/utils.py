import yaml
def read_yaml(file_location):
    yaml_file = {}
    with open(file_location, 'r') as stream:
        try:
            yaml_file = yaml.load(stream, Loader=yaml.FullLoader)
            return yaml_file
        except yaml.YAMLError as exc:
            print(exc)