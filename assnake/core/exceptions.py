class InstanceConfigNotFound(Exception):
    def __init__(self, message="Instance configuration not found or is invalid. Please run 'assnake config init' to configure the instance."):
        super().__init__(message)