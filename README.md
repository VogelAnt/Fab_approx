# Faber approximation
Both scripts use the Faber polynomial approach for the approximation of the matrix exponential 

One relies on the Joukovsky mapping to get the recurrence relation, the other the Schwarz Christoffel mapping 

The matrices we had to look at were usually non hermitian so the whole approach relies on using the symmetry of the spectrum well.  

![spectrum_1](https://user-images.githubusercontent.com/42518184/107871339-c098e380-6ea0-11eb-86d1-73f5fe38429e.png)  

In this case the non hermitian matrix contains eigenvalues whose imaginary part is either 0 or complex conjugated.
So we can exploit the Re-Axis symmetry and use a rectangle based on the numerical range to create a suitable domain,  

![rectangle_spectrum](https://user-images.githubusercontent.com/42518184/107871331-be368980-6ea0-11eb-8179-4a51dbaf7c76.png)  

that is framing the spectrum. With the Faber polynomials in blue :  


![rectangle](https://user-images.githubusercontent.com/42518184/107871330-bd9df300-6ea0-11eb-97d3-19f42b6a8fef.png)
![rectangle_spectrum_1](https://user-images.githubusercontent.com/42518184/107871333-becf2000-6ea0-11eb-8c0a-2f8de5f9ebc8.png)

![rectangle_spectrum](https://user-images.githubusercontent.com/42518184/107871331-be368980-6ea0-11eb-8179-4a51dbaf7c76.png)


![reentrant_faber_10](https://user-images.githubusercontent.com/42518184/107871334-becf2000-6ea0-11eb-9227-e4582a7628df.png)
![reentrant_faber_40](https://user-images.githubusercontent.com/42518184/107871335-bf67b680-6ea0-11eb-8c54-16f22b364d65.png)

An overview of how Schwarz Christoffel mappings work
![SC-mappings](https://user-images.githubusercontent.com/42518184/107871336-c0004d00-6ea0-11eb-92e9-46a97c0dbfd7.png)
![SC-mappings_1](https://user-images.githubusercontent.com/42518184/107871337-c0004d00-6ea0-11eb-8cbb-eca8bd7d19ca.png)

What the extermap type looks like for a random matrix
![extermap](https://user-images.githubusercontent.com/42518184/107871329-bd9df300-6ea0-11eb-81ad-3eb93c873c82.png)



Both approaches preserve the CPTP property of the density matrix as shown below compared to conventional exp(At) solver
![Ergebnis](https://user-images.githubusercontent.com/42518184/107871327-bbd42f80-6ea0-11eb-94b3-1ac963e98c9e.png)
