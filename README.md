
# manifoldGPgrad
## R-package for Bayesian inference on rates on change for spatial processes over compact Riemanian manifolds
Statistical inference on the rates of change of spatial processes has seen some recent developments. Modern spatial data often arise from non-Euclidean domains. This R-package particularly considers spatial processes defined over compact Riemannian manifolds. We develop a computational framework for spatial rates of change for such processes over vector fields. Predictive inference on these rates is devised conditioned on the realized process over the manifold. Manifolds arise as polyhedral meshes in practice. We use a fully Bayesian hierarchical model-based approach to inference on the differential processes for a spatial process over a manifold from partially realized data. 

### Package Structure:
R/: contains required subroutines for spatial modeling and inference with clear documentation for each function and their arguments

inst/extdata: contains the mesh files for the bunny. Load it like so, after installing the package: 

``` r
bunny_path = system.file("extdata", "bun_zipper_res4.ply", package = "manifoldGPgrad")
```

Install the package like so:

``` r
devtools::install_github("arh926/manifoldGPgrad")
```
Below are some examples of what it generates.

## Spatial Process Modeling and inference for gradients on the 2-sphere 
<img width="3228" height="1275" alt="sphere_sim-combined_scale" src="https://github.com/user-attachments/assets/f1c6f6ee-04bf-4e5d-a346-55d6c2555c07" />

## Vector Fields for Directional Derivatives
<img width="2468" height="586" alt="vector_fields_bunny_bl" src="https://github.com/user-attachments/assets/77ffdf64-a4b3-4851-b331-8df1dd078938" />

## Spatial Process Modeling and inference for gradients on the Stanford Bunny
<img width="2710" height="1275" alt="bunny_sim-combined_scale" src="https://github.com/user-attachments/assets/f3d16059-1eef-4ff7-9f9c-62c7a8158152" />

Below is an example working script that helps in reproducing the analysis for the Stanford Bunny. The 2-sphere works similarly.

``` r
############################################
# A plotting function to visualize results #
############################################
# plot_mesh <- function(mesh = NULL,
#                       mesh_obj = NULL,
#                       ncolors = 11,
#                       palette = "Spectral"){
# 
#   interp = as.numeric(interpolate_mesh(mesh = mesh, sample = mesh_obj, method = "spline"))
#   cols = colorRampPalette(RColorBrewer::brewer.pal(ncolors, palette))(100)
#   # scale [0, 1]
#   idx  <- cut((interp - min(interp))/(max(interp) - min(interp)), 100)
# 
#   open3d()
#   shade3d(bunny, color = cols[idx], meshColor = "vertices")
# 
# }
# Note:: Must set options BEFORE loading the rgl package
options(rgl.useNULL = TRUE)
options(rgl.printRglwidget = TRUE)

require(Rvcg)
require(rgl)
require(ggplot2)
require(Matrix)
require(manifoldGPgrad)


# Read the Stanford Bunny
bunny_path = system.file("extdata", "bun_zipper_res4.ply", package = "manifoldGPgrad")
bunny = vcgPlyRead(bunny_path, updateNormals = TRUE)
n = ncol(bunny$vb); n

# Cotangent-Laplacian
L_cot = laplacian_cotangent(mesh = bunny)
L_sym = (L_cot + t(L_cot)) / 2  # just to be safe

eig   = eigen(L_sym, symmetric = TRUE, only.values = FALSE)
ord = order(eig$values) # ascending order: important for the sampling procedure
lambda = eig$values[ord]
phi    = eig$vectors[,ord]


k_eig = 200  # number of eigenpairs
lambda_k = lambda[1:k_eig]
phi_k    = phi[, 1:k_eig]

######################
# Fit a GP to Bunny #
#####################
set.seed(1234)
nsamp = 1e3
coords = sample_points_mesh(mesh = bunny, npoints = nsamp)
tau2 = 1
nu = 2

gp_sample = rnorm(nsamp, 10 * (sin(3 * pi * coords$points[,1]) + sin(3 * pi * coords$points[,2]) + sin(3 * pi * coords$points[,3])), tau2)
gp_sample = gp_sample - mean(gp_sample)


gp_obj = list(points = coords$points,
              values = gp_sample,
              tri_idx = coords$tri_idx,
              bary = coords$bary)

# Not in R-package
plot_mesh(mesh = bunny, mesh_obj = gp_obj)

mc_sp = bayes_sp_manifold(mesh = bunny,
                          lambda = lambda_k,
                          phi = phi_k,
                          y = gp_sample,
                          samp = coords,
                          nu = nu)

zest = spsample(mesh = bunny,
                lambda = lambda_k,
                phi = phi_k,
                y = gp_sample,
                alpha = mc_sp$alpha,
                sig2 = mc_sp$sig2,
                tau2 = mc_sp$tau2,
                samp = coords,
                nu = nu)
mc_sp$z = zest$z
mc_sp$beta0 = zest$beta0

y_mcmc = matrix(NA, nrow = nrow(zest$z), ncol = ncol(zest$z))
for(i in 1:nrow(zest$z)) y_mcmc[i, ] = mc_sp$z[i,] + mc_sp$beta0[i]

sig2_est = round(quantile(mc_sp$sig2, probs = c(0.5, 0.025, 0.975)), 2)
tau2_est = round(quantile(mc_sp$tau2, probs = c(0.5, 0.025, 0.975)), 2)
alpha_est = round(quantile(mc_sp$alpha, probs = c(0.5, 0.025, 0.975)), 2)
z_est = apply(mc_sp$z, 2, median)

yest = z_est + median(mc_sp$beta0)
final_y = cbind(true = gp_sample,
                est = yest,
                lower_ci = apply(mc_sp$z, 2, quantile, probs = 0.025) + quantile(mc_sp$beta0, probs = 0.025),
                upper_ci = apply(mc_sp$z, 2, quantile, probs = 0.975) + quantile(mc_sp$beta0, probs = 0.975))
# y_sd = apply(y_mcmc, 2, sd)

ggplot(final_y, aes(x = est, y = true)) +
  geom_point(color = "black") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "grey70", alpha = 0.5,
              colour = NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(
    x = "Estimated Process",
    y = "True Process",
    title = "",
    subtitle = "",
    caption = ""
  )

est = rbind(sig2_est, tau2_est, alpha_est); est

gp_obj = list(points = coords$points,
              values = yest,
              tri_idx = coords$tri_idx,
              bary = coords$bary)
# Not in R-package
plot_mesh(mesh = bunny, mesh_obj = gp_obj)


#######################
# Gradient Estimation #
#######################
# Grid over bunny
coords_grid = sample_points_mesh(mesh = bunny, npoints = 4e3)
grid = farthest_point_subsample(coords = coords_grid, nkeep = 4e2)

open3d()
shade3d(bunny, color = "grey80", alpha = 0.3)
points3d(grid$points, col = "black", size = 2)

bunny_gradients = gradient_manifold(mesh = bunny, coords = coords,
                                    model = mc_sp, grid = grid,
                                    lambda = lambda_k, phi = phi_k,
                                    d = 2)
grad_est = bunny_gradients$grad_ci[,1]

# Not in R-package
plot_mesh(mesh = bunny, mesh_obj = bunny_gradients$gp_grad_est_obj)

#########
# Truth #
#########
# Truth
# Tangent vector field
mesh = bunny
mesh_norm  = vcgUpdateNormals(mesh)
Vcoords = t(mesh$vb[1:3, ])
F = t(mesh$it)
V_raw = t(mesh_norm$normals)[, 1:3]
V_vert = matrix(0, nrow(V_raw), 3)

ref = c(1, 0, 0)

for(i in seq_len(nrow(Vcoords))){
  n = V_raw[i,]

  v_tan = ref - sum(ref * n) * n
  v_tan = v_tan/sqrt(sum(v_tan^2))

  V_vert[i,] = v_tan
}

u_grid = grid$bary[,1]
v_grid = grid$bary[,2]
w_grid = grid$bary[,3]

v1_grid_idx = F[grid$tri_idx, 1]
v2_grid_idx = F[grid$tri_idx, 2]
v3_grid_idx = F[grid$tri_idx, 3]

# vector field at grid points (Ngrid x 3)
V_grid = u_grid * V_vert[v1_grid_idx, , drop = FALSE] +
  v_grid * V_vert[v2_grid_idx, , drop = FALSE] +
  w_grid * V_vert[v3_grid_idx, , drop = FALSE]

V_grid = V_grid/sqrt(rowSums(V_grid^2))

N_grid = u_grid * V_raw[v1_grid_idx, , drop = FALSE] +
  v_grid * V_raw[v2_grid_idx, , drop = FALSE] +
  w_grid * V_raw[v3_grid_idx, , drop = FALSE]
N_grid = N_grid/sqrt(rowSums(N_grid^2))


# True Gradient
grad_m_grid = cbind(30 * pi * cos(3 * pi * grid$points[,1]),
                     30 * pi * cos(3 * pi * grid$points[,2]),
                     30 * pi * cos(3 * pi * grid$points[,3]))  # Ngrid x 3

dotgN = rowSums(grad_m_grid * N_grid)
grad_m_tan = grad_m_grid - dotgN * N_grid

V_tan = V_grid - rowSums(V_grid * N_grid) * N_grid
V_tan = V_tan/sqrt(rowSums(V_tan^2))

true_grad = rowSums(grad_m_tan * V_grid)

gp_grad_true_obj = list(points = grid$points,
                        values = true_grad,
                        tri_idx = grid$tri_idx,
                        bary = grid$bary)
# Not in R-package
plot_mesh(mesh = bunny, mesh_obj = gp_grad_true_obj)

final = cbind(true_grad, bunny_gradients$grad_ci) # get grad_ci from estimated_bunny_grads_...
sum(apply(final, 1, function(x) ifelse(x[1] >= x[3] & x[1] <= x[4], 1, 0)))/nrow(grid$points) * 100 # 200 eigen pairs:: 93%; 300 eigen pairs:: 95.25%
colnames(final) = c("true" ,"est", "lower_ci", "upper_ci")

ggplot(final, aes(x = est, y = true)) +
  geom_point(color = "black") +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "grey70", alpha = 0.5,
              colour = NA) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14)) +
  labs(
    x = "Estimated Gradient",
    y = "True Gradient",
    title = "",
    subtitle = "",
    caption = ""
  )

```
For questions and bugs, please reach out to: [Email Me](mailto:ah3758@drexel.edu)


