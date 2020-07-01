#!/home/xwu/bin/python
#-*- coding:utf-8 -*-
import argparse
from Sample import *
from multiprocessing import Pool
from Model import *

def feature(args):
	mytime('feature construction start')
	sample_name='NA12878'

	output_file=args.o if args.o!=None else sample_name+'_feature.txt'
	Sample(sample_name,args.b,args.v,args.c,args.o,args.g,args.f)

	mytime('feature construction end')



def train(args):
	targets=np.loadtxt(args.f,skiprows=1,dtype=str,delimiter='\t',)
	del_group=np.array([l for l in targets if l[4]=='DEL'])
	dup_group=np.array([l for l in targets if l[4]=='DUP'])

	#X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
	




	mytime('train end')



def genotype(args):
	mytime('genotype end')

	clf=joblib.load('model.pkl')
	LR.predict(X_test)

	mytime('genotype end')




















parser=argparse.ArgumentParser(prog='sff',usage='python sff.py xxx')
sub_parsers=parser.add_subparsers()

#feature
feature_parser=sub_parsers.add_parser('feature')

feature_parser.add_argument('-b',required=True)
feature_parser.add_argument('-v',required=False)
feature_parser.add_argument('-c',required=True)
feature_parser.add_argument('-o')
feature_parser.add_argument('-g',required=False)
feature_parser.add_argument('-f',required=False)
feature_parser.set_defaults(function=feature)

#train
train_parser=sub_parsers.add_parser('train')
train_parser.add_argument('-f',required=True)
train_parser.add_argument('-m',required=False)
train_parser.add_argument('-o')
train_parser.set_defaults(function=train)


#genotype
genotype_parser=sub_parsers.add_parser('genotype')
genotype_parser.add_argument('-g',required=True)
genotype_parser.add_argument('-b',required=True)
genotype_parser.add_argument('-v',required=True)
genotype_parser.add_argument('-c',required=True)
genotype_parser.add_argument('-o')
genotype_parser.set_defaults(function=genotype)


args=parser.parse_args()
args.function(args)
