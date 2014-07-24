#generate figures of TUD for different clusters, etc
import os
os.chdir('C:\\Users\\Admin\\Documents\\GitHub\\tango\\src')
from tetranucleotideAnalysis import *

# plot all phages
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 phage',True, '..\\figures\\all_phage_TUD.png')
# individual clusters 
#A1
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 A1 phage',True, '..\\figures\\A1_phage_TUD.png', subset=range(0,8))
#A2 
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 A2 phage',True, '..\\figures\\A2_phage_TUD.png', subset=range(8,13))
#B1
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 B1 phage',True, '..\\figures\\B1_phage_TUD.png', subset=range(13,16))
#B2
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 B2 phage',True, '..\\figures\\B2_phage_TUD.png', subset=range(16,18))
#B3
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 B3 phage',True, '..\\figures\\B3_phage_TUD.png', subset=range(18,20))
#B4
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 B4 phage',True, '..\\figures\\B4_phage_TUD.png', subset=range(20,22))
#C1
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 C1 phage',True, '..\\figures\\C1_phage_TUD.png', subset=range(22,28))
#C2
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 C2 phage',True, '..\\figures\\C2_phage_TUD.png', subset=range(28,29))
#D
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 D phage',True, '..\\figures\\D_phage_TUD.png', subset=range(29,35))
#E
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 E phage',True, '..\\figures\\E_phage_TUD.png', subset=range(35,39))
#F1
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 F1 phage',True, '..\\figures\\F1_phage_TUD.png', subset=range(39,47))
#F2
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 F2 phage',True, '..\\figures\\F2_phage_TUD.png', subset=range(47,48))
#G
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 G phage',True, '..\\figures\\G_phage_TUD.png', subset=range(48,50))
#H1
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 H1 phage',True, '..\\figures\\H1_phage_TUD.png', subset=range(50,52))
#H2
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 H2 phage',True, '..\\figures\\H2_phage_TUD.png', subset=range(52,53))
#I
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 I phage',True, '..\\figures\\I_phage_TUD.png', subset=range(53,55))
#Sin
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 Singleton phage',True, '..\\figures\\Singleton_phage_TUD.png', subset=range(55,60))
#representatives from different clusters
plotTUD('hatful60TUD.tsv', 'Hatful et al 2010 representative phage',True, '..\\figures\\rep_phage_TUD.png', subset=[0,8,13,16,18,20,22,28,29,35,39,47,48,50,52,53,55])

