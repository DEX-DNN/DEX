import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

print('Estimate the expected distribution of contributions ... ')
eta=1.5;h=2.5
m = np.linspace(0.001,0.999,1000);
p=1/(np.sqrt(2*3.1415926)*eta*h)*1/(m*(1-m**(1/h)))*np.exp(-(np.log(m**(-1/h)-1)**2/(2*1*eta**2)))

percent_p=np.zeros(20)
percent_e=np.zeros(20)
for it_thresh in range(20):
    percent_p[it_thresh]=(quad(lambda x:1/(np.sqrt(2*3.1415926)*1)*np.exp(-(np.log(x**(-1/h)-1)**2/(2*1*eta**2)))*1/(eta*h*x*(1-x**(1/h))),it_thresh*0.05,1)[0])/(quad(lambda x:1/(np.sqrt(2*3.1415926)*1)*np.exp(-(np.log(x**(-1/h)-1)**2/(2*1*eta**2)))*1/(eta*h*x*(1-x**(1/h))),0,1)[0])
    percent_e[it_thresh]=(quad(lambda x:x*1/(np.sqrt(2*3.1415926)*1)*np.exp(-(np.log(x**(-1/h)-1)**2/(2*1*eta**2)))*1/(eta*h*x*(1-x**(1/h))),it_thresh*0.05,1)[0])/(quad(lambda x:x*1/(np.sqrt(2*3.1415926)*1)*np.exp(-(np.log(x**(-1/h)-1)**2/(2*1*eta**2)))*1/(eta*h*x*(1-x**(1/h))),0,1)[0])
    print('=> %3d percents pixels task over %3d range and %3d percents of contributoins.' % (int(percent_p[it_thresh]*100), int((1-0.05*it_thresh)*100), int(percent_e[it_thresh]*100)))

plt.figure()
ax1=plt.subplot(1,2,1);ax1.set_title('Probability vs. Score');
ax1=ax1.plot(m,p)
ax2=plt.subplot(1,2,2);ax2.set_title('Percent vs. Enegery');
ax2=ax2.plot(percent_p,percent_e)

plt.savefig('Expected_Distributions.jpg')
