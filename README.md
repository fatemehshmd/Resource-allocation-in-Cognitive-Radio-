# Resource-allocation-in-Cognitive-Radio
Resource allocation for underlay DSA Cognitive Radio networks using reinforcement learning (Q-Learning))
Introduction:
-------------
This is a repository for the implementation of cognitive radio network. Since quality measurement of end user plays an ever increasing role in development of the wireless communications toward the 5G era, mean opinion score (MOS) has become a widely used metric, not only because it reflects the subjective quality experience of end users but it also provides a common quality assessment metric for traffic of different types. This paper presents a distributed underlay dynamic spectrum access (DSA) scheme based on MOS which performs integrated traffic management and resource allocation across traffics of dissimilar characteristics (real-time video and data traffic). The presented scheme maximizes the overall MOS through a reinforcement learning for a system where primary users coexist with secondary users accessing the same frequency band of interest, while satisfying a total interference constraint to the primary users. The use of MOS as a common metric allows teaching between nodes carrying different traffic without reducing performance. As a result, the docitive paradigm is applied to the presented scheme to investigate the impact of different docition scenarios on overall MOS where a new comer node being taught by experienced peers with similar and dissimilar traffics. 


Running Scripts:
----------------
The project includes different scripts. The main code is AllTogether_New_Docitive.m, at the beginning of this script, there is a some regarding the simulated scenario (network topology + the number of secondary users and the Q-learning parameters). The action set can be changed in create_state_set.m, and the distance between the users (PU and SU) and their coresponding base stations can be set.

Acknowledgement:
----------------
If you use our model in your research, please cite the following paper:
QoE-Driven Integrated Heterogeneous Traffic Resource Allocation Based on Cooperative Learning for 5G Cognitive Radio Networks,Fatemeh Shah Mohammadi, Andres Kwasinski, IEEE 5G World Forum, 9-11 July 2018, Santa Clara, California, USA.

