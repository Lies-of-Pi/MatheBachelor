# Code of the Delphi model
import torch
import torch.nn as nn
from torch.nn import functional as F
import math


# define GeLU activation function
def gelu(x):
    return 0.5 * x * (1.0 + torch.tanh(math.sqrt(2.0 / math.pi) * (x + 0.044715 * torch.pow(x, 3.0))))

# LayerNorm
class LayerNorm(nn.Module):
    """ LayerNorm but with an optional bias. PyTorch doesn't support simply bias=False """

    def __init__(self, ndim, bias):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(ndim))
        self.bias = nn.Parameter(torch.zeros(ndim)) if bias else None

    def forward(self, input):
        return F.layer_norm(input, self.weight.shape, self.weight, self.bias, 1e-5)

# MLP
class MLP(nn.Module):
    def __init__(self, n_embd):
        super().__init__()
        self.linear1 = nn.Linear(n_embd, 4 * n_embd, bias=True)
        self.linear2 = nn.Linear(4 * n_embd, n_embd, bias=True)

    def forward(self, x):
        x = self.linear1(x)
        x = gelu(x)
        x = self.linear2(x)
        return x

# define Self Attention mechanism by hand
class Attention(nn.Module):
    def __init__(self, n_embd):
        super().__init__()
        self.linear_q = nn.Linear(n_embd, n_embd, bias=False)
        self.linear_k = nn.Linear(n_embd, n_embd, bias=False)
        self.linear_v = nn.Linear(n_embd, n_embd, bias=False)
        self.n_embd = n_embd

    def forward(self, x):
        q = self.linear_q(x)
        k = self.linear_k(x)
        v = self.linear_v(x)
        att = F.softmax(q @ k.transpose(-2, -1) * 1 / math.sqrt(self.n_embd), dim=-1) @ v
        return att       

# internal Block
class Block(nn.Module):
    def __init__(self, n_embd):
        super().__init__()
        self.ln_1 = LayerNorm(n_embd, bias=True)
        self.attn = Attention(n_embd)
        self.ln_2 = LayerNorm(n_embd, bias=True)
        self.mlp = MLP(n_embd)

    def forward(self, x):
        y = self.attn(self.ln_1(x)) 
        x = x + y + self.mlp(self.ln_2(x + y))
        return x

# Next Token Prediction Block
class NextToken(nn.Module):
    def __init__(self):
        super().__init__()
    def forward(self, x):
        return F.softmax(torch.sum(x, dim=0), dim=-1)

# Time-to-Event Prediction Block
class TimeToEvent(nn.Module):
    def __init__(self, T_min=0):
        super().__init__()
        self.T_min = T_min
    def forward(self, x):
        #y = - torch.log(torch.exp(-torch.sum(torch.logsumexp(x,-1))) + self.T_min)
        rates = torch.exp(torch.sum(x, axis=0))
        rate = torch.sum(rates)
        y = torch.min(rates)
        y = torch.max(torch.tensor([y, self.T_min]))
        return y, rate

# Delphi Model
class Delphi(nn.Module):
    def __init__(self, D, n_embd, T_min=0):
        super().__init__()
        self.embedding = nn.Linear(D, n_embd, bias=False)
        self.block = Block(n_embd)
        self.layernorm = LayerNorm(n_embd, bias=True)
        self.linear = nn.Linear(n_embd, D, bias=True)
        self.next_token = NextToken()
        self.time_to_event = TimeToEvent(T_min)
        self.D = D
        self.n_embd = n_embd

    def forward(self, x):
        x = self.embedding(x)
        x = self.block(x)
        x = self.linear(self.layernorm(x))
        y = self.next_token(x)
        t, rate = self.time_to_event(x)
        return y, t, rate