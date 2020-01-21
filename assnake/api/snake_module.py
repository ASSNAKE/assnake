class SnakeModule:
    name = ''
    install_dir = ''
    snakefiles = []
    invocation_commands = []

    def __init__(self, name, install_dir, snakefiles, invocation_commands):
        self.name = name
        self.install_dir = install_dir
        self.snakefiles = snakefiles
        self.invocation_commands = invocation_commands
