import numpy as np

pulls = 1000 #Set this to the desired number of parameter set pulls, we recommend at least 1000 pulls in order to get good statistics
xp = 1  #Projectile: proton=0 , neutron=1
xz = 20 #Nucleus proton number
xa = 40 #Nucleus mass number
xe = 5 #Scattering energy in MeV
numparam = 46 #Number of model parameters, this should stay fixed

#Data arrays for global optical potential parameterizations of each chiral interaction
proton_params=np.array([[ 5.77490e+01,  3.37557e-01,  2.47603e-04,  1.14440e-06,
         1.28964e+01,  1.40872e-01,  2.42242e-04,  1.26886e+00,
         3.10880e-01,  8.78189e-04,  5.50032e-06,  7.79820e-01,
         2.69169e-04,  3.07670e-06,  1.14317e-01,  6.64886e-01,
         4.51031e+00,  3.06178e-01,  5.83818e-04,  1.53636e+01,
         8.98301e-02,  5.47559e-01,  5.61655e+01,  6.73607e-01,
         6.63520e+01,  6.50176e-01,  1.85816e-06,  5.50560e-01,
         2.74471e-01,  1.41153e+01,  4.25837e-01,  2.41938e-03,
         9.44387e-01,  3.72283e-02,  6.86595e-01,  1.31159e-01,
         1.02370e+00,  8.12740e-01,  5.68496e-03,  7.69628e-01,
         1.00135e+01,  8.54929e-03,  1.22918e+00,  8.17408e-01,
         7.74917e-01,  3.17690e-04],
       [ 5.46629e+01,  3.09019e-01,  1.92926e-04,  1.54550e-06,
         1.32845e+01,  1.49706e-01,  2.57806e-04,  1.30974e+00,
         3.90921e-01,  5.09112e-04,  1.98107e-06,  7.71890e-01,
         1.93914e-04,  1.41720e-06,  7.46330e-02,  5.22274e-01,
         3.85483e+00,  2.36678e-01,  4.65075e-04,  1.17927e+01,
         7.67014e-02,  4.76297e-01,  6.07702e+01,  7.75712e-01,
         6.44102e+01,  4.86032e-01,  1.19471e-06,  5.29877e-01,
         3.03935e-01,  1.66306e+01,  3.89183e-01,  1.94573e-03,
         1.30715e+00,  7.71023e-02,  2.16729e+00,  2.87691e-01,
         1.07346e+00,  7.58108e-01,  6.23037e-03,  7.86314e-01,
         9.39118e+00,  8.02939e-03,  1.25967e+00,  8.26821e-01,
         7.80790e-01,  3.74290e-04],
       [ 5.30214e+01,  2.74221e-01,  2.87954e-04,  2.19070e-07,
         1.52556e+01,  1.12644e-01,  1.54952e-04,  1.24417e+00,
         2.10650e-01,  6.23078e-04,  4.69257e-06,  7.44400e-01,
         3.34990e-05,  3.50650e-06,  4.53410e-02,  4.97237e-01,
         2.83618e+00,  1.41413e-01, -1.32877e-07,  9.90006e+00,
         4.80778e-02,  7.52581e-01,  3.31257e+01,  4.59008e-01,
         4.86148e+01,  7.56983e-01,  1.85630e-06,  5.17600e-01,
         2.61587e-01,  1.32263e+01,  3.76659e-01,  1.73911e-03,
         9.95634e-01,  3.04126e-02,  9.47036e-01,  8.88250e-02,
         1.00375e+00,  4.84553e-01,  6.94145e-03,  7.65379e-01,
         1.00957e+01,  8.49747e-03,  1.23265e+00,  7.92527e-01,
         7.69466e-01,  3.71410e-04],
       [ 4.95877e+01,  2.98785e-01, -7.10428e-04,  4.52150e-06,
         1.07123e+01,  1.18676e-01,  4.74972e-04,  1.36490e+00,
         2.58107e-01,  1.04697e-03,  1.18306e-05,  8.01940e-01,
        -8.21600e-05,  1.11010e-05, -1.59040e-02,  5.50789e-01,
         5.13123e+00,  1.95539e-01,  3.83729e-04,  9.26132e+00,
         7.57181e-02,  6.74997e-01,  5.67454e+01,  6.71340e-01,
         7.12548e+01,  5.41396e-01,  7.36243e-07,  5.80005e-01,
         2.98058e-01,  1.78713e+01,  4.09606e-01,  2.47731e-03,
         7.61906e-01,  2.75799e-02,  1.89441e-01,  8.60480e-02,
         9.04490e-01,  0.00000e+00,  0.00000e+00,  8.76141e-01,
         8.86358e+00,  8.46551e-03,  1.34759e+00,  9.75889e-01,
         8.54430e-01,  3.52910e-04],
       [ 5.19440e+01,  2.70587e-01, -7.72760e-04,  3.90610e-06,
         1.22607e+01,  1.40724e-01,  5.75778e-04,  1.34695e+00,
         2.90109e-01,  9.76092e-04,  1.14254e-05,  8.02580e-01,
         1.27713e-04,  9.49250e-06, -1.89850e-02,  5.25134e-01,
         5.60157e+00,  2.33915e-01,  1.96782e-04,  9.73027e+00,
         9.26776e-02,  8.57206e-01,  2.06816e+01,  4.49658e-01,
         3.87414e+01,  7.74680e-01,  1.87261e-06,  6.37886e-01,
         2.14395e-01,  1.00734e+01,  4.92388e-01,  3.64285e-03,
         4.34000e-01,  0.00000e+00,  0.00000e+00,  0.00000e+00,
         7.15000e-01,  0.00000e+00,  0.00000e+00,  6.27500e-01,
         9.59626e+00,  9.03525e-03,  1.33041e+00,  9.53738e-01,
         8.49260e-01,  3.44680e-04]])

neutron_params=np.array([[ 5.71581e+01,  3.28055e-01,  1.93092e-04,  1.27824e-06,
         2.09515e+01,  2.67401e-01,  7.47553e-04,  1.25258e+00,
         3.21132e-01,  6.70967e-04,  4.24950e-06,  7.25804e-01,
         4.16904e-04,  3.86419e-06,  2.14543e-01,  7.58107e-01,
         3.25537e+00,  3.12599e-01,  6.25048e-04,  9.99957e+00,
         7.54228e-02,  5.59196e-01,  6.76906e+01,  6.89906e-01,
         7.64380e+01,  9.33670e-01,  2.57898e-06,  2.19943e-01,
         5.94204e-01,  6.05321e+00,  1.21647e-02,  4.84285e-04,
         2.18277e+00,  4.58225e-02,  3.02741e+00,  1.40682e-01,
         1.26379e+00,  1.48432e+00,  4.20908e-03,  8.20927e-01,
         1.00135e+01,  8.54929e-03,  1.22918e+00,  8.17408e-01,
         7.74917e-01,  3.17690e-04],
       [ 5.40220e+01,  2.99006e-01,  1.41452e-04,  1.66736e-06,
         2.08434e+01,  2.65492e-01,  6.62695e-04,  1.29798e+00,
         3.97229e-01,  5.40576e-04,  1.97630e-06,  7.21033e-01,
         3.94364e-04,  1.74962e-06,  1.70648e-01,  6.24857e-01,
         2.21505e+00,  2.61652e-01,  5.81812e-04,  7.60384e+00,
         4.76306e-02,  4.65605e-01,  8.85857e+01,  8.17832e-01,
         8.52889e+01,  8.32774e-01,  2.55027e-06,  2.04605e-01,
         6.30732e-01,  8.33116e+00,  7.10580e-02,  8.42787e-04,
         1.87864e+00,  4.04378e-02,  3.29653e+00,  1.53675e-01,
         1.30567e+00,  1.47385e+00,  3.91707e-03,  8.33677e-01,
         9.39118e+00,  8.02939e-03,  1.25967e+00,  8.26821e-01,
         7.80790e-01,  3.74290e-04],
       [ 5.22564e+01,  2.62373e-01,  2.05452e-04,  4.39954e-07,
         2.25296e+01,  2.38283e-01,  6.31830e-04,  1.24546e+00,
         2.33753e-01,  6.55036e-04,  5.23030e-06,  7.14386e-01,
         5.33923e-04,  5.24113e-06,  2.98017e-01,  8.17472e-01,
         1.82420e+00,  1.51790e-01,  3.81671e-05,  7.01318e+00,
        -1.01862e-03,  6.08065e-01,  6.73344e+01,  6.28708e-01,
         6.80065e+01,  9.60325e-01,  2.99417e-06,  1.89499e-01,
         6.10633e-01,  1.05308e+01,  8.35152e-02,  7.51665e-04,
         1.03690e+00,  5.08780e-03,  8.87440e-01,  3.52080e-02,
         1.26186e+00,  1.19301e+00,  3.58826e-03,  8.33927e-01,
         1.00957e+01,  8.49747e-03,  1.23265e+00,  7.92527e-01,
         7.69466e-01,  3.71410e-04],
       [ 4.89638e+01,  2.84317e-01, -8.84999e-04,  5.26989e-06,
         1.99247e+01,  3.23092e-01,  1.62662e-03,  1.35913e+00,
         3.07903e-01,  1.31969e-03,  1.57760e-05,  7.71034e-01,
         1.50701e-03,  2.03980e-05,  3.50893e-01,  1.09001e+00,
         3.27876e+00,  2.45623e-01,  7.36418e-04,  8.63975e+00,
        -2.67419e-04,  5.00823e-01,  1.10825e+02,  8.51822e-01,
         1.13970e+02,  7.69755e-01,  4.16265e-06,  2.65264e-01,
         5.66004e-01,  4.64734e+00,  4.74981e-02,  1.09281e-04,
         1.79793e+00,  6.41075e-02,  2.95102e+00,  2.26054e-01,
         1.34526e+00,  1.70187e+00,  2.02361e-03,  8.84760e-01,
         8.86358e+00,  8.46551e-03,  1.34759e+00,  9.75889e-01,
         8.54430e-01,  3.52910e-04],
       [ 5.12707e+01,  2.55931e-01, -9.42339e-04,  4.62002e-06,
         2.09509e+01,  3.30513e-01,  1.65245e-03,  1.34302e+00,
         3.36900e-01,  1.24823e-03,  1.52620e-05,  7.72197e-01,
         1.24094e-03,  1.83160e-05,  3.30972e-01,  1.03029e+00,
         3.70016e+00,  2.80893e-01,  5.25183e-04,  1.00883e+01,
         8.73254e-03,  6.46281e-01,  6.95997e+01,  6.85829e-01,
         9.13809e+01,  1.05856e+00,  5.74254e-06,  3.12911e-01,
         5.33470e-01,  4.41333e+00,  5.41500e-02,  7.73086e-04,
         1.42733e+00,  4.63125e-02,  2.04585e+00,  1.70510e-01,
         1.27478e+00,  1.56846e+00,  1.88844e-03,  8.78288e-01,
         9.59626e+00,  9.03525e-03,  1.33041e+00,  9.53738e-01,
         8.49260e-01,  3.44680e-04]])

#Calculation of mean and covariance matrix for global optical potential parameters
pmean = [np.array(proton_params[:,i]).mean() for i in range(0,numparam)]
nmean = [np.array(neutron_params[:,i]).mean() for i in range(0,numparam)]

pcol , ncol = [] , []
for i in range(numparam):
	pcol.append(proton_params[:,i])
	ncol.append(neutron_params[:,i])
pcol , ncol = np.vstack(pcol) , np.vstack(ncol)
pcov , ncov = np.cov(pcol) , np.cov(ncol)

#Calculation of Woods-Saxon parameters for WLH optical potential
for i in range(pulls):

	px = np.random.multivariate_normal(pmean, pcov)
	nx = np.random.multivariate_normal(nmean, ncov)
	
	if (xp == 0):
		V = px[0] - px[1]*xe + px[2]*xe**2 + px[3]*xe**3 + px[4]*(xa-2*xz)/xa - px[5]*xe*(xa-2*xz)/xa + px[6]*xe**2*(xa-2*xz)/xa
		rv = px[7] - px[8]/(xa**(1/3)) - px[9]*xe + px[10]*xe**2
		av = px[11] - px[12]*xe - px[13]*xe**2 - px[14]*(xa-2*xz)/xa + px[15]*((xa-2*xz)/xa)**2
		W = px[16] + px[17]*xe - px[18]*xe**2 + px[19]*(xa-2*xz)/xa - px[20]*xe*(xa-2*xz)/xa
		rw = px[21] + (px[22] + px[23]*xa)/(px[24] + xa + px[25]*xe) + px[26]*xe**2
		aw = px[27] - (px[28]*xe)/(-px[29] - xe) + px[30]*(xa-2*xz)/xa - px[31]*xe*(xa-2*xz)/xa
		if (xe < 20) and (xa < 100):
			WS = px[32] - px[33]*xe - px[34]*(xa-2*xz)/xa + px[35]*xe*(xa-2*xz)/xa
			rws = px[36] - px[37]/(xa**(1/3)) - px[38]*xe
			aws = px[39]
		else:
			WS = 0.0
			rws = 0.0
			aws = 0.1
		VSO = px[40] - px[41]*xa
		rso = px[42] - px[43]/(xa**(1/3))
		aso = px[44] - px[45]*xa

	if (xp == 1):
		V = nx[0] - nx[1]*xe + nx[2]*xe**2 + nx[3]*xe**3 - nx[4]*(xa-2*xz)/xa + nx[5]*xe*(xa-2*xz)/xa - nx[6]*xe**2*(xa-2*xz)/xa
		rv = nx[7] - nx[8]/(xa**(1/3)) - nx[9]*xe + nx[10]*xe**2
		av = nx[11] + nx[12]*xe - nx[13]*xe**2 - nx[14]*(xa-2*xz)/xa + nx[15]*((xa-2*xz)/xa)**2
		W = nx[16] + nx[17]*xe - nx[18]*xe**2 - nx[19]*(xa-2*xz)/xa - nx[20]*xe*(xa-2*xz)/xa
		rw = nx[21] + (nx[22] + nx[23]*xa)/(nx[24] + xa + nx[25]*xe) + nx[26]*xe**2
		aw = nx[27] - (nx[28]*xe)/(-nx[29] - xe) + nx[30]*(xa-2*xz)/xa - nx[31]*xe*(xa-2*xz)/xa
		if (xe < 40):
			WS = nx[32] - nx[33]*xe - nx[34]*(xa-2*xz)/xa + nx[35]*xe*(xa-2*xz)/xa
			rws = nx[36] - nx[37]/(xa**(1/3)) - nx[38]*xe
			aws = nx[39]
		else:
			WS = 0.0
			rws = 0.0
			aws = 0.1
		if (WS < 0.0): WS = 0.0
		VSO = nx[40] - nx[41]*xa
		rso = nx[42] - nx[43]/(xa**(1/3))
		aso = nx[44] - nx[45]*xa

	paramout = [[V,rv,av,W,rw,aw,WS,rws,aws,VSO,rso,aso]]

	if (i == 0):
		paramouttotal = paramout
	else:
		paramouttotal = np.append(paramouttotal,paramout,axis=0)
np.savetxt('WLH_pn'+ str(xp).zfill(1) +'_z'+ str(xz).zfill(3) +'a'+ str(xa).zfill(3) +'_e'+str(xe).zfill(3)+'_pulls'+str(pulls)+'.txt', np.c_[paramouttotal], fmt='%1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f %1.3f')






#Copyright (c) 2021 Taylor Whitehead
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
