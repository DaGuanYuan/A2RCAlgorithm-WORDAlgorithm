# A2RCAlgorithm-WORDAlgorithm
## Introduction
Pattern synthesis is one of the core problems of modern communication, and it has been widely concerned by people all the time. 

Although the traditional array pattern synthesis method can complete the synthesis of the desired pattern, it lacks certain flexibility, that is, when the target pattern changes, the weight vector of the array manifold must be completely redesigned. Therefore, Here I am trying my best to take flexible and precise array response control algorithm as the research topic, focusing on the newly proposed **_A2RC algorithm_** (an Accurate Array Response Control algorithm for pattern synthesis) and **_WORD algorithm_** (array response control based on Weight vector ORthogonal Decomposition). 

In this program, the corresponding control algorithm has been explored and researched. Also, I am **_the first one_** to use the above algorithm to simultaneously synthesize the main lobe and side lobes of the multi-main lobes target pattern in the uniform linear array and the non-uniform linear array with simultaneous constraints on the main lobe beam ripple and side lobe null steering.

## A2RC Algorithm
See [DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/code_A2RC](https://github.com/DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/tree/master/code_A2RC)

## WORD Algorithm
See [DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/code_WORD](https://github.com/DaGuanYuan/A2RCAlgorithm-WORDAlgorithm/tree/master/code_WORD)

##  Parameters in Pattern Synthesis

### Beamwidth
The **_beamwidth_** of the antenna is the width of the mainlobe. Here we choose 3dB band as beamwidth. 

![image](https://user-images.githubusercontent.com/40145471/129480828-22d59f96-2801-44d3-8af3-2a412ea25491.png)

### Ripple
Mainlobe **_ripple_** usually refers to the fluctuation of the normalized power at the mainlobe within a certain range. Considering that the array pattern is the output of the spatial filter, we define the mainlobe ripple of the array antenna pattern similar to the definition method of the ripple in the filter design.

![image](https://user-images.githubusercontent.com/40145471/129480834-c4632cb3-f9d2-4f4c-a9b0-8d001dbec645.png)

### Sidelobe Level and Nulls
The pattern of all antennas has a maximum value and a minimum value. We call the minimum value(s) the **_null_**, as showed below. In the figure, The pattern is a beampattern of a ULA composed of 6 array elements adjacent to each other with a half-wavelength distance. Considering that in actual work, the scanning range of the mainlobe is between [-60°, 60°]. For this reason, when we adjust the null depth, we only consider the null within the range. In figure below, there are a total of 4 zero nulls, which are located on -40°, -20°, 20°, and 40° respectively. The absolute value that defines the difference between the minimum value of the null and the maximum value of the pattern is called the null depth, as showed below, the null depth at is about -68dB.

![image](https://user-images.githubusercontent.com/40145471/129480837-2e893aed-c370-459a-afa2-99f00a9644ea.png)

However, in fact, it is not convenient to directly control the null depth. Therefore, here we control the null depth by controlling the peak sidelobe levels (PSLs) of the pattern, which is usually called **_Sidelobe levels_**. Specifically, the control of the null depth can be achieved only by adjusting sidelobe levels' normalized power.

## Synthesised Pattern Obtained by A2RC and WORD Algorithm

### Requirements
- ULAs and non-ULAs
- 2 or more Mainlobes
- Ripple < 1dB
- Sidelobe level < -25dB

### Initial Pattern for Uniform Linear Array(ULAs)
Consider a ULA with 18-elements, and the spacing between adjacent elements is half wavelength. We assume that the two main lobes of the beampattern point at -45° and 45° respectively. Here, in order to achieve multiple main lobes, I adopt the following strategy: first make the beampattern point at 45°, and then adjust the normalized power value at -45° to 0dB through A2RC (or WORD algorithm).

![image](https://user-images.githubusercontent.com/40145471/129460549-a58cbb2d-4f64-48a3-97a2-d2c0d492ce26.png)

#### Synthesised Pattern Obtained by A2RC Algorithm for Uniform Linear Array(ULAs)
Considering that the A2RC algorithm multi-point control model stipulates that the mainlobe should be synthesized first, we apply the A2RC multi-point control algorithm to the mainlobe area, as shown in the figure below.

![image](https://user-images.githubusercontent.com/40145471/129460379-4a6d3823-4b9f-4ef5-a8d1-90044a00ca7a.png)

As mentioned, the mainlobe width usually adopts the 3dB bandwidth of the original beam pattern. As shown in the figure above, Δθ<sup>3dB</sup> is about 17.2°, and the main lobe is located in [-48.60°, -41.06°] and [41.06°, 48.60°] (marked in the figure). It can be seen from the figure that the main lobe ripple is controlled within 1dB through repeated iterations of the A2RC algorithm.

Subsequently, we control the sidelobes based on the results of the mainlobe control. Considering that the mainlobe area may change when the side lobes are adjusted, we will continue to adapt both the mainlobes and sidelobes until meeting the final requirements.

![image](https://user-images.githubusercontent.com/40145471/129460385-05388527-ebe7-4c74-a52d-a484ecc13665.png)

After repeated iterations, the target requirements were finally successfully achieved (see figure above). As shown by the dotted line at the ordinate -25dB in the figure, all the sidelobe peak normalized powers have been adjusted to -25dB, meaning that the null depth must be greater than 25dB.

For the main lobe area, considering that the beam pattern of the main lobe area are usually convex upwards, in particular, we calibrated the value of the boundary of the mainlobe area, that is, _L_(dB), the minimum value in the area, which meets the requirements.

|Angle of Arrival(°)|_L_(dB)|
|-------------------|-------|
|-48.60|-0.9566|
|-41.06|-0.9862|
|41.06|-0.8468|
|48.60|-0.9971(_L_<sub>min</sub>)|


#### Synthesised Pattern Obtained by WORD Algorithm for Uniform Linear Array(ULAs)
Firstly, We adapt mainlobe using the WORD multi-point control algorithm.

![image](https://user-images.githubusercontent.com/40145471/129460395-7d40a190-7a3f-4c76-8144-7d72f58aca45.png)

Because of 3dB bandwidth, as showed above, the width of mainlobe area is about 17.2°, and the main lobe is located in [-48.60°, -41.06°] and [41.06°, 48.60°] (Same as A2RC algorithm because of same initial pattern). Mainlobe ripple is adjusted to less than 1dB by repeating the WORD algorithm.

And then, sidelobes would be adapted iteratively until both the mainlobes and sidelobes meet the final requirements.

![image](https://user-images.githubusercontent.com/40145471/129460401-8a460223-8d3e-4ca1-9518-1bf1e2952afc.png)

As we can see, all sidelobes are adjusted to less than -25dB and the boundaries of mainlobes are showed below:

|Angle of Arrival(°)|_L_(dB)|
|-------------------|-------|
|-48.60|-0.9372(_L_<sub>min</sub>)|
|-41.06|-0.3967|
|41.06|-0.5970|
|48.60|-0.9372|
|-45.04|0.5580(_L_<sub>max</sub>)|

It should be noticed that the WORD algorithm allows the normalized power value to exceed 1. Therefore, we should mark both the boundary value (the minimum value in the area) and the maximum value in the area.

### Initial Pattern for Non-Uniform Linear Array
For convenience, we use ULA to design our non-uniform linear array. Consider a ULA with 21 elements, and the spacing between adjacent elements is half a wavelength. We randomly remove the 3<sup>rd</sup>, 7<sup>th</sup>, 9<sup>th</sup>, 13<sup>th</sup>, and 14<sup>th</sup> elements in the ULA to obtain a non-uniform linear array. Through the A2RC algorithm (or WORD algorithm), two main lobes of the nonlinear array beam pattern to -10° and 30° are synthesized respectively, as shown below.

![image](https://user-images.githubusercontent.com/40145471/129484068-32090263-270a-4e2e-a8a5-80cb610f9fa5.png)

#### Synthesised Pattern Obtained by A2RC Algorithm for Non-Uniform Linear Array
Considering that the mainlobe should be synthesized first, here we use the A2RC multi-point control algorithm to the mainlobe area.

![image](https://user-images.githubusercontent.com/40145471/129460518-a2042aee-0d54-4c9a-8d17-f87feda38bac.png)

Since mainlobe width usually adopts the 3dB bandwidth of the original beam pattern, the mainlobe width, which is, Δθ<sup>3dB</sup> is about 7.5°. Also, the main lobe is located between [-12.00°, -8.50°] and [28.36°, 32.34°]. The main lobe ripple is adapted within 1dB through repeated iterations of the A2RC algorithm.

Followed by mainlobe controlling results, we control the sidelobes using WORD algorithm repeatedly until both the mainlobe   and sidelobes meet the final requirements since the mainlobe area may change when the sidelobes are adjusted.

![image](https://user-images.githubusercontent.com/40145471/129460522-b179744d-1c1c-467c-99aa-d9576e598b5b.png)

As showed above, the target requirements were finally successfully achieved. All the sidelobe peak normalized powers have been adjusted to -25dB, which means that the null depth must be greater than 25dB.

Considering that the beampatterns of the mainlobe area are usually convex upwards, we calculate the normalized power of the mainlobe area boundaries, which is, _L_(dB).

|Angle of Arrival(°)|_L_(dB)|
|-------------------|-------|
|-12.00|-0.9396(_L_<sub>min</sub>)|
|-8.50|-0.4335|
|28.36|-0.4474|
|32.34|-0.2829|

All values of _L_ are greater than -1dB, meaning the ripple, which is 0dB - _L_, would less than 1dB.

#### Synthesised Pattern Obtained by WORD Algorithm for Non-Uniform Linear Array
Next, we will use the WORD algorithm to integrate the initial non-uniform linear array multi-mainlobe beampattern into a pattern that meets the ## Requirements. We first use the WORD algorithm to control the main lobe area of the beam pattern (the result is shown in the figure below). The main lobe area and 3dB bandwidth of the beam pattern are the same as the results of the A2RC algorithm, and it can be seen intuitively that the main lobe ripple is synthesized below 1dB.

![image](https://user-images.githubusercontent.com/40145471/129460527-9a5ade33-fcca-44ef-90d9-8107ca6ad38c.png)

After the mainlobe control is completed, we start to control the sidelobe area. Similarly, we need to constantly adjust the mainlobe and sidelobes to counter the changes that occur in other areas outside the adjustment area. As shown in the figure below, after continuous adjustment, all the nulls have been systhesized to less than -25dB.

![image](https://user-images.githubusercontent.com/40145471/129460532-acf24a79-1989-40d2-8a26-206d0013ace8.png)

Considering that the WORD algorithm allows the normalized power value to exceed 1, therefore, for the result of the mainlobe area, we also need to consider the part above 0dB. In particular, we mark both the boundary value (the minimum value in the area) and the maximum value in the area,

|Angle of Arrival(°)|_L_(dB)|
|-------------------|-------|
|-12.00|-0.5776|
|-8.50|-0.9967(_L_<sub>min</sub>)|
|28.36|-0.3171|
|32.34|-0.5372|
|29.84|-0.0280(_L_<sub>max</sub>)|

## Future Works
Although this project has done some exploration and research on the precise beam emission control algorithm of the array antenna pattern, and made some meagre contributions, there are still many issues worthy of in-depth discussion.

- First of all, the multi-point control algorithm of the A2RC algorithm and the WORD algorithm are realized by iterating the single-point control algorithm. But in fact, in the iterative process of synthesizing **_BOTH_** the main lobe and the side lobes at the same time, it may happen that when the side lobes are synthesized, the main lobe changes and does not meet the index requirements; when the main lobe is synthesized, the side lobes may change The lobe will also change and not meet the index requirements, so the algorithm will enter an infinite loop that **_cannot be converged_**. The author believes that this is the inevitable drawback of obtaining the multi-point control model through single-point control iteration, because the energy of the received signal of the array is conserved. When one point is adjusted, it will inevitably affect other points. Therefore, the author believes that **_the multi-point control model should be directly modeled theoretically in future research to avoid such an iterative infinite loop._**

- In addition, the selection of key parameters and is not perfect in the A2RC algorithm and the WORD algorithm. The A2RC algorithm gives two criteria for the selection of parameters, namely formula (44) and formula (52) in paper A2RC. However, it is worth noting that formula (44) is based on global search, which will inevitably **_lead to a huge amount of calculation_**. Although formula (52) reduces the amount of calculation, it may lead to **_Severe distortion_** of the beam pattern. Similarly, the WORD algorithm is not very rigorous in the selection of parameters. The selection criterion (29) in paper WORD is very **_empirical_**, that is, the criterion cannot guarantee the establishment of formula (44). Hence, there is still a lot of mathematical work that can be done in the parameter selection criteria of the two algorithms.

- Finally, the constraint on the normalized power in the WORD algorithm does not strictly follow the definition. That is to say, the normalized power is allowed to **_be greater than 1_** in the WORD algorithm, which is caused by the particularity of its iterative algorithm (not updated θ<sub>0</sub>). What is interesting is that if the normalized power is subject to certain constraints and limited to the range of 0 to 1, the convergence of the WORD algorithm will often be **_very poor_**. Therefore, how to **_update the iterative process variable_** of the WORD algorithm is still worthy of discussion and research.
