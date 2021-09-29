import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from matplotlib.patches import Rectangle

def compute_proper_time(time,alpha):
    # Compute proper time
    proper_time = np.zeros(len(time))
    for n in range(1,len(time)):
        proper_time[n] = proper_time[n-1] + 0.5*(time[n]-time[n-1])*(alpha[n]+alpha[n-1])

    return proper_time

def draw_rectangle_and_connecting_lines(ax1_in,ax2_in,rec_left,rec_right,rec_bot,rec_top,connect):

    # Draw rectangle
    con = Rectangle((rec_left,rec_bot),rec_right-rec_left,rec_top-rec_bot,edgecolor='black',facecolor='none',lw=0.5,ls='--')
    ax1_in.add_artist(con)

    # Draw connecting lines from plot 1 to plot 2
    if connect == 'bottom_to_top':
        xys_1 = [(rec_left,rec_top),(rec_right,rec_top)]
        xys_2 = [(rec_left,rec_bot),(rec_right,rec_bot)]
        for i in range(2):
            con = ConnectionPatch(xyA=xys_1[i], xyB=xys_2[i], coordsA="data", coordsB="data",
                                  axesA=ax2_in, axesB=ax1_in,color="black",ls='--',lw=0.5)
            ax1_in.add_artist(con)
    elif connect == 'left_to_right':
        xys_1 = [(rec_left,rec_top),(rec_left,rec_bot)]
        xys_2 = [(rec_right,rec_top),(rec_right,rec_bot)]
        for i in range(2):
            con = ConnectionPatch(xyA=xys_1[i], xyB=xys_2[i], coordsA="data", coordsB="data",
                                  axesA=ax1_in, axesB=ax2_in,color="black",ls='--',lw=0.5)
            ax1_in.add_artist(con)

# Output file name
outfile = "lapse_self_similarity.png"

# Load weak data
t_w,alp_w,sf_w = np.loadtxt("out_weak.dat").T

# Load strong data
t_s,alp_s,sf_s = np.loadtxt("out_strong.dat").T

fig = plt.figure()

linewidth = 1
color_w   = 'blue'
color_s   = 'orange'

# Top panel
ax1 = plt.subplot(2, 1, 1)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w,label=r"$\eta_{\rm weak} = 0.3364266156435$")
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--',label=r"$\eta_{\rm strong} = 0.3364266156436$")
plt.legend(loc=1,markerfirst=False)
plt.xlim(-0.4,6.4)
plt.ylim(-0.05,0.75)
plt.xticks([0,1,2,3,4,5,6],['0','1','2','3','4','5','6'])
plt.yticks([0,0.2,0.4,0.6],['0.0','0.2','0.4','0.6'])

# Bottom right panel
xl2 = [+5.06,+5.34]
yl2 = [-0.03,+0.23]
ax2 = plt.subplot(2, 2, 4)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w)
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--')
plt.xlim(xl2[0],xl2[1])
plt.ylim(yl2[0],yl2[1])
plt.xticks([5.08,5.2,5.32],['5.08','5.20','5.32'])
plt.yticks([0,0.1,0.2],['0.0','0.1','0.2'])

# Bottom left panel
xl3 = [+5.3,+5.32]
yl3 = [-0.02,+0.12]
ax3 = plt.subplot(2, 2, 3)
plt.grid(ls=':')
plt.plot(t_w,alp_w,lw=linewidth,c=color_w)
plt.plot(t_s,alp_s,lw=linewidth,c=color_s,ls='--')
plt.xlim(xl3[0],xl3[1])
plt.ylim(yl3[0],yl3[1])
plt.xticks([5.302,5.31,5.318],['5.302','5.310','5.318'])
plt.yticks([0,0.06,0.12],['0.00','0.06','0.12'])

draw_rectangle_and_connecting_lines(ax1,ax2,xl2[0],xl2[1],yl2[0],yl2[1],'bottom_to_top')
draw_rectangle_and_connecting_lines(ax2,ax3,xl3[0],xl3[1],yl3[0],yl3[1],'left_to_right')

# Set labels
fig.text(0.5, 0.03, r"$t$", ha='center', va='center')
fig.text(0.03, 0.5, r"$\alpha_{\rm central}$", ha='center', va='center', rotation='vertical')

plt.savefig(outfile,dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)

# Compute proper time
tauw = compute_proper_time(t_w,alp_w)

# Compute logarithmic proper time
tau_star = 1.22958674
N        = len(tauw[tauw<=tau_star])
xiw      = np.zeros(N)
sfxiw    = np.zeros(N)
for i in range(N):
    xiw[i]   = -np.log(np.abs(tau_star-tauw[i]))
    sfxiw[i] = sf_w[i]

# Begin plotting
fig,ax = plt.subplots(ncols=3,figsize=(15,5),dpi=300,sharey=True)

label_fontsize = 22
ticks_fontsize = 20

for i in range(3):
    ax[i].grid()

# Plot 1: Central scalar field as a function of coordinate time
ax[0].set_ylabel(r"$\varphi_{\rm central}$",fontsize=label_fontsize)
ax[0].set_xlabel(r"$t$",fontsize=label_fontsize)
ax[0].plot(t_w,sf_w,'blue')
ax[0].set_xticks([0,1,2,3,4,5])
ax[0].set_xticklabels(["0","1","2","3","4","5"],fontsize=ticks_fontsize)
ax[0].set_yticks([-0.6,-0.4,-0.2,0,0.2,0.4,0.6])
ax[0].set_yticklabels(["-0.6","-0.4","-0.2","0.0","0.2","0.4","0.6"],fontsize=ticks_fontsize)

# Plot 2: Central scalar field as a function of proper time
ax[1].set_xlabel(r"$\tau$",fontsize=label_fontsize)
ax[1].plot(tauw,sf_w,'blue')
ax[1].set_xticks([0,0.3,0.6,0.9,1.2])
ax[1].set_xticklabels(["0.0","0.3","0.6","0.9","1.2"],fontsize=ticks_fontsize)

# Plot 3: Central scalar field as a function of proper time
ax[2].set_xlabel(r"$\xi \equiv - \ln\left|\tau_{*}-\tau\right|$",fontsize=label_fontsize)
ax[2].plot(xiw,sfxiw,'blue')
ax[2].set_xticks([0,3,6,9,12,15])
ax[2].set_xticklabels(["0","3","6","9","12","15"],fontsize=ticks_fontsize)

outfile = "central_scalar_field.png"
plt.tight_layout()
plt.savefig(outfile,dpi=300,bbox_inches='tight',facecolor='white')
plt.close(fig)
