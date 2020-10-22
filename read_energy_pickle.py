import pickle, os

def pickle2text():

    pfile = open('weig.dat','rb')
    data = pickle.load(pfile)
    pfile.close()
    
    fname = open('weig.txt','w')
    print(data,file=fname)
    fname.close()

    return
    
def main():
    pickle2text()
    return
    
main()
