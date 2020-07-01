#!/home/xwu/bin/python
#-*- coding:utf-8 -*-

from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
import numpy as np
import os

class Model():
	def __init__(self,x,y,m=None,pkl='model.pkl'):
		if m!=None:
			self._m=m
		elif os.path.exists(pkl):
			self._m=joblib.load(pkl)
			LR.predict(X_test)
		elif x!=None and y!=None:
			self._m=self.Train(x,y)
		else:
			print 'model error!'
			exit(1)

	def Train(self,x,y,pkl='model.pkl',coefs_file='coefs.txt'):
		LR = LogisticRegression(C=100000.0, penalty='l2', tol=0.0001,solver='liblinear')
		LR.fit(x,y)
		print 'The training score is',LR.score(x,y)
		joblib.dump(LR,pkl)
		scale=abs(LR.intercept_)
		coefs=LR.coef_[0]/scale
		intercept=LR.intercept_/scale
		f=open(coefs_file,'w')
		f.write('coefs:\n')
		for i in range(len(coefs)-1):
			f.write(str(coefs[i])+'\t')
		f.write(str(coefs[-1])+'\n')
		f.write('intercept:\n'+str(intercept)+'\n')
		f.close()
		return LR
	
	def Test(self,x):
		if self._m!=None:
			return self._m.predict(x)
		else:
			print 'model is missing!'











