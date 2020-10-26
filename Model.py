#!/home/xwu/bin/python
#-*- coding:utf-8 -*-

from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
import numpy as np
import os

class Model():
	def __init__(self,filename,suffix='model.pkl'):
		data=self.loadFeature(filename)
		if data[0]=='train':
			self.Train(data,suffix,feature_name)
		else:
			self.Test(data,suffix)

		if os.path.exists(pkl):
			self._m=joblib.load(pkl)
			LR.predict(X_test)
		elif x!=None and y!=None:
			self._m=self.Train(x,y)
		else:
			print 'model error!'
			exit(1)

	def Train(self,data,suffix='model.pkl',feature_name,coefs_file='coefs.txt'):
		del_model = LogisticRegression(C=100000.0, penalty='l2', tol=0.0001,solver='liblinear')
		dup_model = LogisticRegression(C=100000.0, penalty='l2', tol=0.0001,solver='liblinear')
		del_model.fit(data[1][0],data[1][1])
		dup_model.fit(data[2][0],data[2][1])
		print 'The training score of deletion group is',del_model.score(x,y)
		print 'The training score of duplication group is',dup_model.score(x,y)
		joblib.dump(del_model,'del_'+suffix)
		joblib.dump(dup_model,'dup_'+suffix)
		#normalization
		model_dict={'del_model':del_model,'dup_model':dup_model}
		f=open(coefs_file,'w')
		f.write('#\t'+'\t'.join(feature_name)+'\tintercept\n')
		for m in model_dict:
			f.write(m+'\t')
			scale=abs(model_dict[m].intercept_)
			coefs=del_model.coef_[0]/scale
			intercept=LR.intercept_/scale
			for i in range(len(coefs)):
				f.write(str(coefs[i])+'\t')
			f.write(str(intercept)+'\n')
		f.close()

	
	def Test(self,x):
		if self._m!=None:
			return self._m.predict(x)
		else:
			print 'model is missing!'

	def loadFeature(filename):
		data=np.loadtxt(filename,skiprows=0,dtype=str,delimiter='\t',)
		col_name=data[0]
		targets=data[1:]
		del_group=np.array([l for l in targets if l[4]=='DEL'])
		dup_group=np.array([l for l in targets if l[4]=='DUP'])
		if 'category' in col_name:
			del_y=del_group[:,6]
			del_x=del_group[:,7:]
			dup_y=dup_group[:,6]
			dup_x=dup_group[:,7:]
			return 'train',(del_x,del_y),(dup_x,dup_y),col_name[8:]
		else:
			del_x=del_group[:,5:]
			dup_x=dup_group[:,5:]
			return 'test',del_x,dup_x










