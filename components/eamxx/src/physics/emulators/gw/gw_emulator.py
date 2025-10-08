from model_definition import ANN_CNN

model = None

def init (nlevs,num_hidden_layers,dict_filename):
    global model

    model = ANN_CNN(idim=nlevs,odim=nlevs,hdim=num_hidden_layers,stencil=1,dropout=0)
    model.load_state_dict(dict_filename)
    model.eval()

def run (dt,u,v,w):

    pass
