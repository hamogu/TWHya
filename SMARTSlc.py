import atpy
tab = atpy.Table('/data/hguenther/TWHya/SMARTSlc.dat', type = 'ascii')

filters = ['B','V','R','I','J','H']

fig = plt.figure(figsize=(4,3))
ax = fig.add_subplot(111)
for i, filt in enumerate(filters):
    ax.errorbar(tab['JD']-2455650, tab[filt]+i, tab[filt+'_err'], label=filt)

ax.set_xlabel('JD - 2455650')
ax.set_ylabel('mag [arb. zero point]')
ax.legend()
ax.set_xlim([0,80])
ax.set_ylim([6,-1])
fig.subplots_adjust(top = .96, bottom = .16, left=.15, right=.95)
plt.draw()
fig.savefig('/data/hguenther/Dropbox/my_articles/TW_Hya_UV/SMARTSlc.png', bb_inches='tight')

