# Stand Alone - prototyping for label encoding and comparison
# Aldo Saavedra - Health Hackathon Oct 2016
#module to incode each label into binary value
#ATCGN
# A 0001  
# C 0010
# G 0100
# T 1000
# N 1111

labels = {1:'A',2:'C',4:'G',8:'T',15:'N'}

def EncodeSequence(sequence):
    #encode each label of the sequence into a binary value

    enc_sequence=[]
    for label in sequence:

        enc_label = 0
        if (label == 'A'):
            enc_label = 1
        elif (label == 'C'):
            enc_label = 2
        elif (label == 'G'):
            enc_label = 4
        elif (label == 'N'):
            enc_label = 8
        else:
            enc_label = 15  # T

        enc_sequence.append(enc_label)
    
    return enc_sequence

def Compare_Two_Sequence(seq_under_test,std_seq):
    # will comapare the two sequences after converting the labels
    # to numbers ..
    # returns the position where the labels are different

    enc_seq_ut = EncodeSequence(seq_under_test)
    enc_std_seq = EncodeSequence(std_seq)

    if (len(enc_seq_ut) != len(enc_std_seq)):
        print 'Error - the length of the sequences do not match'
        return 

    pos={}
    j=0;
    for i in range(len(enc_seq_ut)):
        if enc_seq_ut[i] != enc_std_seq[i]:
            pos[j]=[i,labels[enc_std_seq[i]]]
            j = j + 1

    return pos



std_sequence='GTCAATGTCACTGAGAGAACTCAGAGAGGAATTATGAGAAAAGCAAACCGGAGTATTTTCAATAGCAAAATTCTGTAACTTTTCATCAGTTGCTGCCCCTCTGTCTGGTATGTCTTTGGATGACTGGGGAAAAGTGGATTGTTTCTGAAGTATGGGTTTAGGCTGACCTCGATTTATTGGCTGCTTTGCAATAGCTTGTGTCTTATTAGCTGATTGTTGGTTGGAGGTTAGTTCTGGT'

test_sequence='GTCTATGTCACTGAGAGAACTCAGAGAGGAATTATGAGAAAAGCAAACCGGAGTATTTTCAATAGCAAAATTCTGTAACTTTTCATCAGTTGCTGCCCCTCTGTCTGGTATGTCTTTGGATGACTGGNNAAAAGTGGATTGTTTCTGAAGTATGGGTTTAGGCTGACCTCGATTTATTGGCTGCTTTGCAATAGCTTGTGTCTTATTAGCTGATTGTTGGTTGGAGGTTAGTTCTGGT'



enc_seq = Compare_Two_Sequence(test_sequence,std_sequence)

print enc_seq



