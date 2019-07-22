#!/home/xwu/bin/python
#-*- coding:utf-8 -*-
import argparse
from Sample import *
from multiprocessing import Pool
from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
import numpy as np

def feature(args):
	mytime('feature construction start')
	sample_name='NA12878'

	output_file=args.o if args.o!=None else sample_name+'_feature.txt'
	Sample(sample_name,args.b,args.v,args.c,args.o,args.g,args.f)

	mytime('feature construction end')



def train(args):
	mytime('train start')
	f=open('NA12878_feature.txt','r')
	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)
	LR = LogisticRegression(C=1.0, penalty='l1', tol=0.01,solver='liblinear',multi_class='auto')
	LR.fit(X_train,y_train)
	print LR.score(X_test,y_test)
	joblib.dump(LR,'model.pkl')




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
