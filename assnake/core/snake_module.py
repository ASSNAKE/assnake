class SnakeModule:
    name = ''
    install_dir = ''
    snakefiles = []
    invocation_commands = []
    initialization_commands = []
    wc_configs = []
    initialization_commands = []
    results = []

    def __init__(self, name, install_dir, snakefiles, invocation_commands, initialization_commands = [], wc_configs = [], results = []):

        self.name = name
        self.install_dir = install_dir
        self.snakefiles = snakefiles
        self.invocation_commands = invocation_commands
        self.initialization_commands = initialization_commands
        self.wc_configs = wc_configs
        self.results = results