class Extractor():
    def __init__(self, name):
        super().__init__()
        self.activation = {}
        self.name = name
        
    def get_activation(self):
        name = self.name
        def hook(model, input, output):
            self.activation[name] = output.detach()
        return hook
            