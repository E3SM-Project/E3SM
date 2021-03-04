import matplotlib
import matplotlib.pyplot as plt
x = [16,32,64,96,192]
y0 = [264.,161.5,124.2,94.3,83.]
y  = [i/60.0 for i in y0]
fig, ax = plt.subplots()
ax.scatter(x,y)
ax.set_title("Performance Test Chart")
ax.set_xlabel('Number of processors')
ax.set_ylabel('Wall clock time [min]')
xlable = [str(i) for i in x]
#ax.set_xticklabels(xlabels)
plt.xticks(x,xlable)
fig.savefig('performance_test.png')


