# A2RCAlgorithm-WORDAlgorithm
## Introduction
Pattern synthesis is one of the core problems of modern communication, and it has been widely concerned by people all the time. 

Although the traditional array pattern synthesis method can complete the synthesis of the desired pattern, it lacks certain flexibility, that is, when the target pattern changes, the weight vector of the array manifold must be completely redesigned. Therefore, Here I am trying my best to take flexible and precise array response control algorithm as the research topic, focusing on the newly proposed **_A2RC algorithm_** (an Accurate Array Response Control algorithm for pattern synthesis) and **_WORD algorithm_** (array response control based on Weight vector ORthogonal Decomposition). 

In this program, the corresponding control algorithm has been explored and researched. Also, I am **_the first one_** to use the above algorithm to simultaneously synthesize the main lobe and side lobes of the multi-main lobes target pattern in the uniform linear array and the non-uniform linear array with simultaneous constraints on the main lobe beam ripple and side lobe null steering.

## A2RC Algorithm
See [DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/code_A2RC](https://github.com/DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/tree/master/code_A2RC)

## WORD Algorithm
See [DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/code_WORD](https://github.com/DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/tree/master/code_WORD)

## Synthesised Pattern Obtained by A2RC and WORD Algorithm

### Initial Pattern for Uniform Linear Array(ULAs)

![image](https://user-images.githubusercontent.com/40145471/129460549-a58cbb2d-4f64-48a3-97a2-d2c0d492ce26.png)

#### A2RC Algorithm for Uniform Linear Array(ULAs)
Mainlobe control:

![image](https://user-images.githubusercontent.com/40145471/129460379-4a6d3823-4b9f-4ef5-a8d1-90044a00ca7a.png)

Both mainlobe and sidelobe control:

![image](https://user-images.githubusercontent.com/40145471/129460385-05388527-ebe7-4c74-a52d-a484ecc13665.png)


#### WORD Algorithm for Uniform Linear Array(ULAs)
Mainlobe control:

![image](https://user-images.githubusercontent.com/40145471/129460395-7d40a190-7a3f-4c76-8144-7d72f58aca45.png)

Both mainlobe and sidelobe control:

![image](https://user-images.githubusercontent.com/40145471/129460401-8a460223-8d3e-4ca1-9518-1bf1e2952afc.png)

### Initial Pattern for Non-Uniform Linear Array

#### A2RC Algorithm for Non-Uniform Linear Array
Mainlobe control:

![image](https://user-images.githubusercontent.com/40145471/129460518-a2042aee-0d54-4c9a-8d17-f87feda38bac.png)

Both mainlobe and sidelobe control:

![image](https://user-images.githubusercontent.com/40145471/129460522-b179744d-1c1c-467c-99aa-d9576e598b5b.png)

#### WORD Algorithm for Non-Uniform Linear Array(ULAs)
Mainlobe control:

![image](https://user-images.githubusercontent.com/40145471/129460527-9a5ade33-fcca-44ef-90d9-8107ca6ad38c.png)

Both mainlobe and sidelobe control:

![image](https://user-images.githubusercontent.com/40145471/129460532-acf24a79-1989-40d2-8a26-206d0013ace8.png)

## Conclusion

## Future Works
Although this project has done some exploration and research on the precise beam emission control algorithm of the array antenna pattern, and made some meagre contributions, there are still many issues worthy of in-depth discussion.

- First of all, the multi-point control algorithm of the A2RC algorithm and the WORD algorithm are realized by iterating the single-point control algorithm. But in fact, in the iterative process of synthesizing **_BOTH_** the main lobe and the side lobes at the same time, it may happen that when the side lobes are synthesized, the main lobe changes and does not meet the index requirements; when the main lobe is synthesized, the side lobes may change The lobe will also change and not meet the index requirements, so the algorithm will enter an infinite loop that **_cannot be converged_**. The author believes that this is the inevitable drawback of obtaining the multi-point control model through single-point control iteration, because the energy of the received signal of the array is conserved. When one point is adjusted, it will inevitably affect other points. Therefore, the author believes that **the multi-point control model should be directly modeled theoretically in future research to avoid such an iterative infinite loop.**

- In addition, in fact, the selection of key parameters and is not perfect in the A2RC algorithm and the WORD algorithm. The A2RC algorithm gives two criteria for the selection of parameters, namely formula (44) and formula (3-49) in the paper of A2RC. However, it is worth noting that formula (44) is based on global search, which will inevitably lead to a huge amount of calculation. Although formula (3-49) reduces the amount of calculation, it may lead to Severe distortion of the beam pattern. Similarly, the WORD algorithm is not very rigorous in the selection of parameters. The selection criterion (3-87) is very empirical, that is, the criterion cannot guarantee the establishment of formula (44). Therefore, the author feels that there is still a lot of mathematical work that can be done in the parameter selection criteria of the two algorithms.

