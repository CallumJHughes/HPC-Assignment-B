import matplotlib.pyplot as plt

per = [33.6, 15.0, 9.8, 1.0, 0.1]

sub = ['[3]', '[4]', '[5]', '[6]', '[7]']

#sub = [1.0,1.5,1.6,1.7,1.8,1.9,2.0]

#per = [417.56483328935633,417.56497749795926,417.56499546879490,417.56501138504018,417.56502555200746,417.56503827800100,530.04653173320878]

plt.bar(sub,per)

#plt.ylim(0,100)

plt.ylabel("Percentage time/%")

plt.xlabel("Subprocess")

plt.title("Time spent in each subprocess of program")

plt.show()

