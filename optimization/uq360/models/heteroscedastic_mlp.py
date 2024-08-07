import torch
import torch.nn.functional as F
from uq360.models.noise_models.heteroscedastic_noise_models import GaussianNoise

class GaussianNoiseMLPNet(torch.nn.Module):

    def __init__(self, num_features, num_outputs, num_hidden_1, num_hidden_2, num_hidden_3):
        super(GaussianNoiseMLPNet, self).__init__()
        self.fc_1 = torch.nn.Linear(num_features, num_hidden_1)
        self.fc_2 = torch.nn.Linear(num_hidden_1, num_hidden_2)
        self.fc_3 = torch.nn.Linear(num_hidden_2, num_hidden_3)
        self.fc_mu = torch.nn.Linear(num_hidden_3, num_outputs)
        self.fc_log_var = torch.nn.Linear(num_hidden_3, num_outputs)
        self.noise_layer = GaussianNoise()

    def forward(self, x):
        x = F.relu(self.fc_1(x))
        x = F.relu(self.fc_2(x))
        x = F.relu(self.fc_3(x))
        mu = self.fc_mu(x)
        log_var = self.fc_log_var(x)
        return mu, log_var

    def loss(self, y_true=None, mu_pred=None, log_var_pred=None):
        return self.noise_layer.loss(y_true, mu_pred, log_var_pred, reduce_mean=True)
    