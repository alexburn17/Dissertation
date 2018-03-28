#The evolution of reproductive restraint through social communication - Review

**P. Alexander Burnham**

March 6, 2018

**Paper Summary:**In this paper, the authors tackle the age old question in evolutionary biology related to the evolution of altruism. More precisely, they examined the evolutionary history of communication related to reproduction and reproductive restraint using agent-based modeling techniques. The authors show that cooperation through communication resulted in favorable outcomes that led to a stable population. The authors make make three predictions: **1)** reproductive restraint, when communicated, is a favorable adaptation in the face of overcrowding, **2)**  Individuals that communicate individuals that do not communicate in a population and **3)** once that population is communicating, it is not susceptible to invasion to non-communicative individuals. The modeling technique used by the authors was to create a cellular automaton where there were two organism in an antagonistic relationship (i.e predator-prey, parasite-host etc.) termed "consumers" and "hosts" respectively. Hosts can reproduce into neighboring cells at a certain probability g, be  confronted with a consumed with a probability, tau, and by a consumer with a probability, v. When the host population overcrowds (individual surrounded on all sides by other hosts), local signals were sent out (communication) to which could tell neighbors to decrease reproduction by some value, delta. Both delta and tau were subject to mutation. Natural, without communication, both values reach a moderate equilibrium.

**Experiment 1:** Is communication of reproductive restraint a favorable adaptation in the face of overcrowding? Delta became more negative than a random model (0) in communication based models meaning that evolution favors communication to reduce reproduction during overcrowding. 

**Experiment 2:** Are individuals that communicate more likely to persist than individuals who don't. In a population of individuals not able to adapt to overcrowding, individuals who could adapt were introduced. When introduced, 18 out of 20 invasion was successful. When non-communicating individuals were introduced under the same parameters, they were not as successful in their invasion.

**Experiment 3:** This is the reverse of experiment 2. Can non communicating individual invade a communicating population? Not one invasion in 141,369 attempts was successful. 

**Critique:**This paper did an excellent job in explaining yet another useful implementation of cellular automata in modeling systems in fields with open questions. Using this modeling technique to examine the advantages of altruistic behavior, a contentious arena in evolutionary biology, is an interesting one, which takes some of the complexity out of the problem and reduces it to see simple parts. The structure of the model and the rules described clearly painted a picture and showed the power of communication in a population.

**Questions and observations:**

* What is a specially extended environment in terms of this paper and the models therein and why is it an assumption for there three predictions? 
* In experiment 2, thy talk about introductions and runs. Runs are only mentioned for the treatment runs (18/20) but not for the control. What is the difference?
* Not a question but it looks like Charles Goodnight, an evolutionary biologist in our biology department was acknowledged in this paper!