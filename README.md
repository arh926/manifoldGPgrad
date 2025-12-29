
# manifoldGPgrad
## R-package for Bayesian inference on rates on change for spatial processes over compact Riemanian manifolds
Statistical inference on the associated rates of change has seen some recent developments. Modern spatial data often arise from non-Euclidean domains. This R-package particularly considers spatial processes defined over compact Riemannian manifolds. We develop a computational framework for spatial rates of change for such processes over vector fields. Predictive inference on these rates is devised conditioned on the realized process over the manifold. Manifolds arise as polyhedral meshes in practice. We use a fully Bayesian hierarchical model-based approach to inference on the differential processes for a spatial process over a manifold from partially realized data. 

Install the package like so:

```{r}
devtools::install_github("arh926/manifoldGPgrad")
```
Below are some examples of what it generates.

<img width="3228" height="1275" alt="sphere_sim-combined_scale" src="https://github.com/user-attachments/assets/f1c6f6ee-04bf-4e5d-a346-55d6c2555c07" />

<img width="2468" height="586" alt="vector_fields_bunny_bl" src="https://github.com/user-attachments/assets/77ffdf64-a4b3-4851-b331-8df1dd078938" />

<img width="2710" height="1275" alt="bunny_sim-combined_scale" src="https://github.com/user-attachments/assets/f3d16059-1eef-4ff7-9f9c-62c7a8158152" />



