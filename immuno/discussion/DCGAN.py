import torch
import torch.nn as nn
from torch.utils.data import Dataset,DataLoader
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import Normalizer, StandardScaler,MinMaxScaler,RobustScaler
from sklearn.metrics import roc_curve, auc, precision_recall_curve, confusion_matrix, f1_score,accuracy_score
import re
import collections




class Discriminator_peptide(nn.Module):
    def __init__(self,channels=1,features=16,embed=12):
        super(Discriminator_peptide,self).__init__()
        self.net = nn.Sequential(
            # N * channel * 10 * embed
            nn.Conv2d(channels,features,kernel_size=(3,embed)),
            nn.LeakyReLU(0.2),
            # N * features * 8 * 1
            nn.Conv2d(features,features*2,kernel_size=(5,1)),
            nn.BatchNorm2d(features*2),
            nn.LeakyReLU(0.2),
            # N * features2 * 4 * 1
            nn.Conv2d(features*2,features*4,kernel_size=(3,1)),
            nn.BatchNorm2d(features*4),
            nn.LeakyReLU(0.2),
            # N * features4 * 2 * 1
            nn.Conv2d(features*4,features*8,kernel_size=(2,1)),
            # N * features8 *1 * 1
            nn.BatchNorm2d(features*8),
            nn.LeakyReLU(0.2),
        )

    def forward(self,x):
        return self.net(x)


class Discriminator_MHC(nn.Module):
    def __init__(self,channels=1,features=16,embed=12):
        super(Discriminator_MHC,self).__init__()
        self.net = nn.Sequential(
            # N * channel * 46 * embed
            nn.Conv2d(channels,features,kernel_size=(15,embed)),
            nn.LeakyReLU(0.2),
            # N * features * 32 * 1
            nn.Conv2d(features,features*2,kernel_size=(9,1),dilation=(2,1)),
            nn.BatchNorm2d(features*2),
            nn.LeakyReLU(0.2),
            # N * features2 * 16 * 1
            nn.Conv2d(features*2,features*4,kernel_size=(5,1),dilation=(2,1)),
            nn.BatchNorm2d(features*4),
            nn.LeakyReLU(0.2),
            # N * features4 * 8 * 1
            nn.Conv2d(features*4,features*8,kernel_size=(3,1),dilation=(2,1)),
            nn.BatchNorm2d(features*8),
            nn.LeakyReLU(0.2),
            # N * features8 * 4 * 1
            nn.Conv2d(features*8,features*8,kernel_size=(4,1)),
            nn.BatchNorm2d(features*8),
            nn.LeakyReLU(0.2),
            # N * features8 * 1 * 1
        )

    def forward(self,x):
        return self.net(x)


class Discriminator(nn.Module):
    def __init__(self,d_peptide,d_MHC):
        super(Discriminator,self).__init__()
        self.d_peptide = Discriminator_peptide(channels=1,features=16,embed=12)
        self.d_MHC = Discriminator_MHC(channels=1,features=16,embed=12)
        self.fc1 = nn.Linear(256,128)
        self.fc2 = nn.Linear(128,1)
        self.relu = nn.ReLU()
        self.sigmoid = nn.Sigmoid()

    def forward(self,x1,x2):   # x1 is peptide input(N,1,10,12), x2 is MHC input(N,1,46,12)
        out1 = self.d_peptide(x1)   # N,feature8,1,1    # N,128,1,1
        out2 = self.d_MHC(x2)    # N, feature8,1,1      # N,128,1,1
        out1 = out1.reshape(out1.shape[0],-1)     # N,128
        out2 = out2.reshape(out2.shape[0],-1)     # N, 128
        out = torch.cat([out1,out2],dim=1)   # N, 256
        out = self.relu(self.fc1(out))
        out = self.fc2(out)
        out = self.sigmoid(out)
        return out


class Generator_peptide(nn.Module):
    def __init__(self,channels=1,features=16,embed=12):
        super(Generator_peptide,self).__init__()
        self.net = nn.Sequential(
            # N * 64 * 1 * 12
            nn.ConvTranspose2d(64,features*2,kernel_size=(4,1)),
            nn.BatchNorm2d(features*2),
            nn.ReLU(),
            # N * 32 * 4 * 12
            nn.ConvTranspose2d(features*2,features,kernel_size=(4,1),stride=(2,1),padding=(1,0)),
            nn.BatchNorm2d(features),
            nn.ReLU(),
            # N * 16 * 8 * 12
            nn.ConvTranspose2d(features,channels,kernel_size=(3,1)),
            nn.Tanh(),
            # N * 1 * 10 * 12
        )

    def forward(self,x):
        return self.net(x)

class Generator_MHC(nn.Module):
    def __init__(self,channels=1,features=16,embed=12):
        super(Generator_MHC,self).__init__()
        self.net = nn.Sequential(
            # N * 256 * 1 * 12
            nn.ConvTranspose2d(256,features*8,kernel_size=(4,1)),
            nn.BatchNorm2d(features*8),
            nn.ReLU(),
            # N * 128 * 4 * 12
            nn.ConvTranspose2d(features*8,features*4,kernel_size=(4,1),stride=(2,1),padding=(1,0)),
            nn.BatchNorm2d(features*4),
            nn.ReLU(),
            # N * 64 * 8 * 12
            nn.ConvTranspose2d(features*4,features*2,kernel_size=(4,1),stride=(2,1),padding=(1,0)),
            nn.BatchNorm2d(features*2),
            nn.ReLU(),
            # N * 32 * 16 * 12
            nn.ConvTranspose2d(features*2,features,kernel_size=(4,1),stride=(2,1),padding=(1,0)),
            nn.BatchNorm2d(features),
            nn.ReLU(),
            # N * 16 * 32 * 12
            nn.ConvTranspose2d(features,1,kernel_size=(8,1),dilation=(2,1)),
            nn.Tanh()
            # N * 1 * 46 * 12
        )

    def forward(self,x):
        return self.net(x)

class Generator(nn.Module):
    def __init__(self,g_peptide,g_MHC):
        super(Generator,self).__init__()
        self.g_peptide = Generator_peptide()
        self.g_MHC = Generator_MHC()
    def forward(self,x1,x2):
        fake_p = self.g_peptide(x1)   # N,64,1,12
        fake_M = self.g_MHC(x2)    # N,256,1,12
        return fake_p,fake_M   # fake_p will be N,1,10,12, fake_M will be N,1,46,12, the same as Discrimator input


class mydataset(Dataset):
    def __init__(self,input1,input2):
        self.input1 = input1
        self.input2 = input2

    def __len__(self):
        return len(self.input1)

    def __getitem__(self,idx):
        x1 = torch.from_numpy(self.input1[idx])
        x2 = torch.from_numpy(self.input2[idx])
        return (x1,x2)

# now preproces some data
def aaindex(peptide,after_pca):

    amino = 'ARNDCQEGHILKMFPSTWYV-'
    matrix = np.transpose(after_pca)   # [12,21]
    encoded = np.empty([len(peptide), 12])  # (seq_len,12)
    for i in range(len(peptide)):
        query = peptide[i]
        if query == 'X': query = '-'
        query = query.upper()
        encoded[i, :] = matrix[:, amino.index(query)]

    return encoded


def peptide_data_aaindex(peptide,after_pca):   # return numpy array [10,21,1]
    length = len(peptide)
    if length == 10:
        encode = aaindex(peptide,after_pca)
    elif length == 9:
        peptide = peptide[:5] + '-' + peptide[5:]
        encode = aaindex(peptide,after_pca)
    encode = encode.reshape(-1,encode.shape[0], encode.shape[1])
    return encode


def hla_data_aaindex(hla, dic_inventory, hla_type,after_pca):    # return numpy array [46,21,1]
    dic = paratope_dic(hla)
    try:
        seq = dic[hla_type]
    except KeyError:
        hla_type = rescue_unknown_hla(hla_type, dic_inventory)
        seq = dic[hla_type]
    encode = aaindex(seq,after_pca)
    encode = encode.reshape(-1,encode.shape[0], encode.shape[1])
    return encode


def construct_aaindex(ori, hla, dic_inventory,after_pca):
    series = []
    for i in range(ori.shape[0]):
        peptide = ori['peptide'].iloc[i]
        hla_type = ori['HLA'].iloc[i]
        immuno = np.array(ori['immunogenecity'].iloc[i]).reshape(1,-1)   # [1,1]

        encode_pep = peptide_data_aaindex(peptide,after_pca)    # [10,12]

        encode_hla = hla_data_aaindex(hla, dic_inventory, hla_type,after_pca)   # [46,12]
        series.append((encode_pep, encode_hla, immuno))
    return series


def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {'A': dicA, 'B': dicB, 'C': dicC}

    for hla in inventory:
        type_ = hla[4]  # A,B,C
        first2 = hla[6:8]  # 01
        last2 = hla[8:]  # 01
        try:
            dic[type_][first2].append(last2)
        except KeyError:
            dic[type_][first2] = []
            dic[type_][first2].append(last2)

    return dic


def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    #print(hla)
    if not big_category.get(first2) == None:
        small_category = big_category.get(first2)
        distance = [abs(int(last2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(first2) + str(optimal)
    else:
        small_category = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small_category]
        optimal = min(zip(small_category, distance), key=lambda x: x[1])[0]
        return 'HLA-' + str(type_) + '*' + str(optimal) + str(big_category[optimal][0])

def paratope_dic(hla):
    df = hla
    dic = {}
    for i in range(df.shape[0]):
        hla = df['hla'].iloc[i]
        paratope = df['paratope'].iloc[i]
        dic[hla] = paratope
    return dic

def add_X(array):
    me = np.mean(array)
    array = np.append(array, me)
    return array


def read_index(path):
    with open(path, 'r') as f:
        data = f.readlines()
        array = []
        for line in data:
            line = line.lstrip(' ').rstrip('\n')
            line = re.sub(' +', ' ', line)

            items = line.split(' ')
            items = [float(i) for i in items]
            array.extend(items)
        array = np.array(array)
        array = add_X(array)
        Index = collections.namedtuple('Index',
                                       ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
                                        'T', 'W', 'Y', 'V', 'X'])
        I = Index._make(array)
    return I, array  # namedtuple


def read_all_indices():
    table = np.empty([21, 566])
    for i in range(566):
        if len(str(i)) == 1:
            ii = '00' + str(i)
        elif len(str(i)) == 2:
            ii = '0' + str(i)
        else:
            ii = str(i)

        NA_list_str = ['472', '473', '474', '475', '476', '477', '478', '479', '480', '481', '520', '523', '524']
        NA_list_int = [int(i) for i in NA_list_str]
        if ii in NA_list_str: continue

        path = '/Users/ligk2e/Desktop/NeoAntigenWorkflow/immunogenecity/AAindex1/index{0}.txt'.format(ii)

        _, array = read_index(path)

        table[:, i] = array
    table = np.delete(table, NA_list_int, 1)
    return table


def scaling(table):  # scale the features
    table_scaled = RobustScaler().fit_transform(table)
    return table_scaled


def wrapper_read_scaling():
    table = read_all_indices()
    table_scaled = scaling(table)
    return table_scaled


def pca_get_components(result):
    pca= PCA()
    pca.fit(result)
    result = pca.explained_variance_ratio_
    sum_ = 0
    for index,var in enumerate(result):
        sum_ += var
        if sum_ > 0.95:
            return index    # 25 components



def pca_apply_reduction(result):
    pca = PCA(n_components=12)  # or strictly speaking ,should be 26, since python is 0-index
    new = pca.fit_transform(result)
    return new



def pull_peptide_aaindex(dataset):
    result = np.empty([len(dataset),1,10,12])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][0]
    return result


def pull_hla_aaindex(dataset):
    result = np.empty([len(dataset),1,46,12])
    for i in range(len(dataset)):
        result[i,:,:,:] = dataset[i][1]
    return result

def pull_label_aaindex(dataset):
    result = np.empty([len(dataset), 1])
    for i in range(len(dataset)):
        result[i, :] = dataset[i][2]
    return result






if __name__ == '__main__':
    # real data
    ori = pd.read_csv('../immuno/data/shuffle_training_test.txt',sep='\t')
    positive = ori.loc[ori['immunogenecity']==1]
    positive.set_index(pd.Index(np.arange(positive.shape[0])),inplace=True)

    # pre-processing
    table_scaled = wrapper_read_scaling()   # [21,553]
    after_pca = pca_apply_reduction(table_scaled)   # [21,12]
    hla = pd.read_csv('hla2paratopeTable_aligned.txt',sep='\t',header=None,names=['hla','paratope'])
    inventory = hla['hla']
    dic_inventory = dict_inventory(inventory)
    dataset = construct_aaindex(positive, hla, dic_inventory,after_pca)   # [ (1,10,12),(1,46,12),(1,1)   ]
    input1 = pull_peptide_aaindex(dataset)
    input2 = pull_hla_aaindex(dataset)
    label_ = pull_label_aaindex(dataset)




    # training
    lr = 0.0005
    batch_size = 128
    num_epochs = 30
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    dataset_ = mydataset(input1,input2)
    dataloader = DataLoader(dataset_,batch_size=batch_size,shuffle=True)

    netD1 = Discriminator_peptide()
    netD2 = Discriminator_MHC()
    netD = Discriminator(netD1,netD2)

    netG1 = Generator_peptide()
    netG2 = Generator_MHC()
    netG = Generator(netG1,netG2)

    optimizerD = torch.optim.Adam(netD.parameters(),lr=lr,betas=(0.5,0.999))
    optimizerG = torch.optim.Adam(netG.parameters(),lr=lr,betas=(0.5,0.999))

    netG.train()
    netD.train()

    criterion = nn.BCELoss()



    # start to train
    for epoch in range(num_epochs):
        lossD_list = []
        lossG_list = []
        for (x1,x2) in dataloader:
            x1 = x1.float().to(device)
            x2 = x2.float().to(device)

            # noise should be re-sample every mini-batch
            noise1 = torch.randn(batch_size, 64, 1, 12).to(device)
            noise2 = torch.randn(batch_size, 256, 1, 12).to(device)

            ## training Discriminator : max log(D(x)) + log(1-D(G(z)))

            # first give all real instances
            netD.zero_grad()
            label = (torch.ones(len(x1)) * 0.9).to(device)
            output = netD(x1,x2).reshape(-1)    # from N,1 to N
            lossD_real = criterion(output,label)

            # second, generate all fake instances to train discriminator
            fake = netG(noise1,noise2)
            label = (torch.ones(batch_size) * 0.1).to(device)
            output = netD(fake[0].detach(),fake[1].detach()).reshape(-1)
            lossD_fake = criterion(output,label)

            lossD = lossD_real + lossD_fake
            lossD.backward()
            optimizerD.step()
            lossD_list.append(lossD.item())

            ## train Generator: max log(D(G(z)))
            netG.zero_grad()
            label = torch.ones(batch_size).to(device)
            output = netD(fake[0],fake[1]).reshape(-1)
            lossG = criterion(output,label)
            lossG.backward()
            optimizerG.step()
            lossG_list.append(lossG.item())

        lossD_mean = sum(lossD_list)/len(lossD_list)
        lossG_mean = sum(lossG_list)/len(lossG_list)
        print('Epoch{0} -- Loss D: {1:.4f} -- Loss G: {2:.4f}'.format(epoch+1,lossD_mean,lossG_mean))


    # save the model and initial noise1 and noise2
    torch.save(netG1.state_dict(),'GAN/netG1.pth')
    torch.save(netG2.state_dict(),'GAN/netG2.pth')
    torch.save(netG.state_dict(),'GAN/netG.pth')



    # generate pseudo-positive samples

    pseudo_p_array = []
    pseudo_MHC_array = []
    for i in range(100):
        noise1 = torch.randn(batch_size, 64, 1, 12).to(device)
        noise2 = torch.randn(batch_size, 256, 1, 12).to(device)
        pseudo_p,pseudo_MHC = netG(noise1,noise2)

        pseudo_p = pseudo_p.detach().cpu().numpy()
        pseudo_MHC = pseudo_MHC.detach().cpu().numpy()
        pseudo_p_array.append(pseudo_p)
        pseudo_MHC_array.append(pseudo_MHC)

    pseudo_p_total = np.concatenate(pseudo_p_array,axis=0)
    pseudo_MHC_total = np.concatenate(pseudo_MHC_array,axis=0)

    import pickle
    with open('/Users/ligk2e/Desktop/tmp/pseudo_p_total.p','wb') as f1:
        pickle.dump(pseudo_p_total,f1)
    with open('/Users/ligk2e/Desktop/tmp/pseudo_MHC_total.p','wb') as f2:
        pickle.dump(pseudo_MHC_total,f2)
















