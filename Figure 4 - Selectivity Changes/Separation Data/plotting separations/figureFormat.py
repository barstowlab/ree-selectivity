fig = plt.figure(figsize=(3.25,2.5))
ax = fig.gca()

plt.setp(ax.get_xticklabels(), fontsize=5)
plt.setp(ax.get_yticklabels(), fontsize=5)
ax.tick_params(axis='both', length=3)
plt.grid(linestyle='--')

leg = ax.legend(prop={"size":5})
ax.set_title('$\Delta$SO_0625 and $\Delta$SO_3385 Increase Yb Biosorption', fontsize=7, fontweight='bold')
ax.set_ylabel(rees[k] + ' Biosorption (nmols)', fontsize=7, fontweight='bold')
ax.set_xlabel('La Biosorption (nmols)', fontsize=7, fontweight='bold')  
fig.tight_layout(pad=0.1)


plt.savefig('clean_del_seps.svg')
plt.savefig('clean_del_seps.eps')