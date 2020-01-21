class SnakeModule:
    name = ''
    install_dir = ''
    snakefiles = []
    invocation_commands = []
    wc_configs = []

    def __init__(self, name, install_dir, snakefiles, invocation_commands, wc_configs=[]):
        self.name = name
        self.install_dir = install_dir
        self.snakefiles = snakefiles
        self.invocation_commands = invocation_commands
        self.wc_configs = wc_configs
