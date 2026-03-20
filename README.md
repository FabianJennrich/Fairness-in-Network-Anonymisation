# Fairness-in-Network-Anonymisation
Implementations of the anonymization algorithms used in the master's thesis on fairness in network anonymisation. This code is written and tested on python 3.11. Works by using PipelineAnonEval.ipynb or by running PipelineAnonEval.py and giving the location of the graph to be anonymized when prompted. The interactANONET.py file is unused in the pipeline, and contains the python functions written for interacting with outputs from ANONET code.

The following Python packages are required: pandas, networkx, datetime, os, random, numpy, itertools, scipy, microagg1d, collections, sklearn.

## Citations
The citations for the papers from which the anonymization methods are implemented are below:

ANONET: de Jong, Rachel G., Mark P. J. van der Loo, and Frank W. Takes. "The anonymization problem in social networks." arXiv preprint arXiv:2409.16163 (2024). doi: https://doi.org/10.48550/arXiv.2409.16163

ANNM: S. Hamzehzadeh and S. M. Mazinani, ‘‘ANNM: A new method for adding noise nodes which are used recently in anonymization methods in social networks,’’ Wireless Pers. Commun., vol. 107, no. 4, pp. 1995–2017, Aug. 2019.

DPFT: H. Zhu, X. Zuo, and M. Xie, ‘‘DP-FT: A differential privacy graph gen- eration with field theory for social network data release,’’ IEEE Access, vol. 7, pp. 164304–164319, 2019.

edgeEnt: J. Yan, L. Zhang, Y. Tian, G. Wen, and J. Hu, ‘‘An uncertain graph approach for preserving privacy in social networks based on important nodes,’’ in Proc. Int. Conf. Netw. Netw. Appl. (NaNA), Oct. 2018, pp. 107–111.

ENDP: K. R. Macwan and S. J. Patel, ‘‘Node differential privacy in social graph degree publishing,’’ Procedia Comput. Sci., vol. 143, pp. 786–793, Jan. 2018.

FiM: F. Beato, M. Conti, and B. Preneel, “Friend in the middle (FiM): Tackling de-anonymization in social networks,” in Proc. 5th Int. Workshop Security Soc. Netw., San Diego, CA, USA, 2013, pp. 279–284.

kDegSeq: J. Casas-Roma, J. Herrera-Joancomartí, and V. Torra, ‘‘K-degree anonymity and edge selection: Improving data utility in large networks,’’ Knowl. Inf. Syst., vol. 50, no. 2, pp. 447–474, Feb. 2017.
microagg1d for the UMGA: Felix I. Stamm and Michael T Schaub. Faster optimal univariate microaggregation. Transactions on Machine Learning Research, 2024.

KDVEM: T.Ma,Y.Zhang,J.Cao,J.Shen,M.Tang,Y.Tian,A.Al-Dhelaan,and M. Al-Rodhaan, ‘‘KDVEM: A k-degree anonymity with vertex and edge modification algorithm,’’ Computing, vol. 97, no. 12, pp. 1165–1184, 2015.

linkMirage: C. Liu and P. Mittal, “Linkmirage: Enabling privacy-preserving ana- lytics on social relationships,” in Proc. NDSS, San Diego, CA, USA, 2016, pp. 1–15.

linkPerturb: P. Mittal, C. Papamanthou, and D. Song, “Preserving link privacy in social network based systems,” in Proc. 20th Annu. Netw. Distrib. Syst. Security Symp. (NDSS), San Diego, CA, USA, 2013, pp. 1–15.

OLP: T. Wu, G. Ming, X. Xian, W. Wang, S. Qiao, and G. Xu, ‘‘Structural predictability optimization against inference attacks in data publishing,’’ IEEE Access, vol. 7, pp. 92119–92136, 2019.

PBCN: H. Huang, D. Zhang, F. Xiao, K. Wang, J. Gu, and R. Wang, ‘‘Privacy- preserving approach PBCN in social network with differential pri- vacy,’’ IEEE Trans. Netw. Service Manage., vol. 17, no. 2, pp. 931–945, Jun. 2020.

pseudoNode: Chester S et al (2013) Why Waldo befriended the dummy? k-Anonymization of social networks with pseudo-nodes. Soc Netw Anal Min 3(3):381–399

randomWalk: Y. Guo, Z. Liu, Y. Zeng, R. Wang, and J. Ma, ‘‘Preserving privacy for hubs and links in social networks,’’ in Proc. Int. Conf. Netw. Netw. Appl. (NaNA), Oct. 2018, pp. 263–269.

upperApprox: S. Kumar and P. Kumar, ’’Upper approximation based privacy preserving in online social networks,’’ Expert Syst. Appl., vol. 88, pp. 276–289, Dec. 2017.

subgraphDP: B. P. Nguyen, H. Ngo, J. Kim, and J. Kim, ‘‘Publishing graph data with subgraph differential privacy,’’ in Proc. Int. Workshop Inf. Secur. Appl. Jeju-do, South Korea: Springer, 2015, pp. 134–145.

twoStep: K. Liu and E. Terzi, “Towards identity anonymization on graphs,” in Proc. ACM SIGMOD Int. Conf. Manag. Data, Vancouver, BC, Canada, 2008, pp. 93–106.	
