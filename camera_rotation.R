# Pinhole camera rotation homography
# www.overfitting.net
# https://www.overfitting.net/

library(Rcpp)
library(tiff)


# Auxiliar functions (the projection transformation equations are implemented in cpp files)

bresenham_line <- function(x0,y0,x1,y1){
    x0<-as.integer(x0);y0<-as.integer(y0);x1<-as.integer(x1);y1<-as.integer(y1)
    dx<-abs(x1-x0); sx<-if(x0<x1)1L else -1L
    dy<--abs(y1-y0); sy<-if(y0<y1)1L else -1L
    err<-dx+dy
    xs<-integer(); ys<-integer()
    repeat{
        xs<-c(xs,x0); ys<-c(ys,y0)
        if(x0==x1 && y0==y1) break
        e2<-2*err
        if(e2>=dy){ err<-err+dy; x0<-x0+sx}
        if(e2<=dx){ err<-err+dx; y0<-y0+sy}
    }
    cbind(x=xs,y=ys)
}

draw_bresenham <- function(img, x0, y0, x1, y1, linewidth = 1, colour = c(1,1,0)) {
    M <- dim(img)[1]
    N <- dim(img)[2]
    pts <- bresenham_line(x0, y0, x1, y1)
    
    radius <- (linewidth - 1) / 2
    rmax <- ceiling(radius)
    
    for(i in seq_len(nrow(pts))) {
        x <- pts[i,1]
        y <- pts[i,2]
        
        for(dx in -rmax:rmax) {
            for(dy in -rmax:rmax) {
                # Disco
                if(dx*dx + dy*dy <= radius*radius) {
                    xx <- x + dx
                    yy <- y + dy
                    
                    if(xx >= 1 && xx <= N &&
                       yy >= 1 && yy <= M) {
                        img[yy, xx, 1] <- colour[1]
                        img[yy, xx, 2] <- colour[2]
                        img[yy, xx, 3] <- colour[3]
                    }
                }
            }
        }
    }
    
    img
}

# Add grid to input image
add_grid <- function(img, n_gridx = 12, linewidth = 4, colour = c(1, 1, 0), diagonal = TRUE) {
    # For 3:2 format images, n_gridx = 12 -> 12x8 squares grid
    # For 4:3 format images, n_gridx = 16 -> 16x12 squares grid
    
    stopifnot(length(dim(img)) == 3, dim(img)[3] == 3)
    M <- dim(img)[1]
    N <- dim(img)[2]
    step <- N / n_gridx
    
    # Si hay un número par de celdas, una línea pasa exactamente por el centro
    has_center_line <- (n_gridx %% 2 == 0)
    
    make_lines <- function(L, step, center_line = TRUE) {
        
        centre <- (L + 1) / 2
        max_steps <- ceiling((L / 2) / step)
        
        if (center_line) {
            multipliers <- seq(-max_steps, max_steps, by = 1)
        } else {
            multipliers <- seq(-max_steps, max_steps, by = 1) + 0.5
        }
        
        pos <- centre + multipliers * step
        pos <- pos[pos >= 1 & pos <= L]
        sort(unique(round(pos)))
    }
    
    expand <- function(v, max_val) {
        expanded <- unlist(lapply(v, function(x) {
            x - floor((linewidth - 1) / 2) + 0:(linewidth - 1)
        }))
        unique(sort(expanded[expanded >= 1 & expanded <= max_val]))
    }
    
    cols <- make_lines(N, step, center_line = has_center_line)
    rows <- make_lines(M, step, center_line = has_center_line)

    cols <- expand(cols, N)
    rows <- expand(rows, M)
    
    # Líneas verticales
    img[, cols, 1] <- colour[1]
    img[, cols, 2] <- colour[2]
    img[, cols, 3] <- colour[3]
    
    # Líneas horizontales
    img[rows, , 1] <- colour[1]
    img[rows, , 2] <- colour[2]
    img[rows, , 3] <- colour[3]
    
    # Diagonales
    if (diagonal) {
        img <- draw_bresenham(img, x0 = 1, y0 = 1, x1 = N, y1 = M, linewidth = linewidth, colour = colour)
        img <- draw_bresenham(img, x0 = N, y0 = 1, x1 = 1, y1 = M, linewidth = linewidth, colour = colour)
    }
    
    img
}


# Calculo de la matriz de rotación para rectificar un trapezoide a rectángulo
# 
# H: alto de la imagen en píxeles
# W: ancho de la imagen en píxeles
# fl_FF_mm: distancia focal en equivalente FF (mm)
# pts: Matriz o data.frame 4x2 con las esquinas (X, Y) en píxeles
#      Orden: 1=Top-Left, 2=Bottom-Left, 3=Bottom-Right, 4=Top-Right
#      Criterio: el píxel superior izquierdo de la imagen se indexa como (1,1)
#                creciendo la coordenada X hacia la derecha y la coordenada Y hacia abajo
get_rectifying_rotation <- function(H, W, fl_FF_mm, pts) {
    
    # 1. Parámetros de cámara
    diag_mm <- sqrt(36.0^2 + 24.0^2)
    diag_pixel <- sqrt(H^2 + W^2)
    f_pixel <- fl_FF_mm * (diag_pixel / diag_mm)
    
    # Centro óptico adaptado para indexación en base 1
    cx <- (W + 1.0) / 2.0
    cy <- (H + 1.0) / 2.0
    
    K_inv <- matrix(c(
        1/f_pixel, 0,         -cx/f_pixel,
        0,         1/f_pixel, -cy/f_pixel,
        0,         0,         1
    ), nrow=3, byrow=TRUE)
    
    # 2. Funciones auxiliares
    cross_prod <- function(a, b) {
        c(a[2]*b[3] - a[3]*b[2],
          a[3]*b[1] - a[1]*b[3],
          a[1]*b[2] - a[2]*b[1])
    }
    normalize <- function(v) { v / sqrt(sum(v^2)) }
    
    # Convertir puntos a coordenadas homogéneas
    P <- cbind(as.matrix(pts), 1)
    
    # 3. Calcular puntos de fuga (ORDEN: TL, BL, BR, TR)
    # Horizontal (eje X): Intersección de la línea superior (TL-TR) y la inferior (BL-BR)
    l_top <- cross_prod(P[1,], P[4,])
    l_bot <- cross_prod(P[2,], P[3,])
    v1_homo <- cross_prod(l_top, l_bot) 
    
    # Vertical (eje Y): Intersección de la línea izquierda (TL-BL) y la derecha (TR-BR)
    l_left <- cross_prod(P[1,], P[2,])
    l_right <- cross_prod(P[4,], P[3,])
    v2_homo <- cross_prod(l_left, l_right) 
    
    # 4. De-proyectar a rayos 3D
    d1 <- normalize(as.vector(K_inv %*% v1_homo))
    d2 <- normalize(as.vector(K_inv %*% v2_homo))
    
    # 5. Corrección de signos (orientación de cámara)
    if (d1[1] < 0) d1 <- -d1
    if (d2[2] < 0) d2 <- -d2
    
    # 6. Forzar ortogonalidad (Gram-Schmidt)
    d2_ortho <- d2 - sum(d1 * d2) * d1
    d2_ortho <- normalize(d2_ortho)
    
    # 7. Eje Z (Normal al plano)
    d3 <- cross_prod(d1, d2_ortho)
    d3 <- normalize(d3)
    
    # 8. Matriz de rotación final
    return(matrix(c(d1, d2_ortho, d3), nrow=3, byrow=TRUE))
}



#######################################
# 0. COMPILAR FUNCIÓN HOMOGRAFÍA EN C++

sourceCpp("camera_rotation.cpp")



#######################################
# 1. ROTACIÓN BÁSICA

# Cargar una imagen de prueba (devuelve un array H x W x 3 con valores de 0 a 1)
img <- readTIFF("street12mm.tif") 
img=add_grid(img, n_gridx = 12, colour = c(1,1,0), linewidth = 4)
writeTIFF(img, "street12mm_grid.tif")


# Crear Matriz de Rotación

deg2rad <- function(deg) {
    return(deg * pi / 180)
}

# 1. Yaw/Paneo: rotación sobre el eje Y
# Equivale a girar la cámara hacia la izquierda o la derecha (como mirar a los lados en Google Street View)
# Modifica las coordenadas X y Z de los rayos de luz, pero mantiene intacta la altura Y

theta_y <- deg2rad(10)  # 10º

R_yaw <- matrix(c(
    cos(theta_y), 0, sin(theta_y),
    0, 1,         0,
   -sin(theta_y), 0, cos(theta_y)
), nrow = 3, byrow = TRUE)


# 2. Pitch/Cabeceo: rotación sobre el eje X
# Equivale a inclinar la cámara hacia arriba o hacia abajo
# Afecta a las coordenadas Y y Z, manteniendo el centro horizontal fijo

theta_x <- deg2rad(5)  # 5º

R_pitch <- matrix(c(
    1,  0,             0,
    0,  cos(theta_x), -sin(theta_x),
    0,  sin(theta_x),  cos(theta_x)
), nrow = 3, byrow = TRUE)


# 3. Roll/Alabeo: rotación sobre el eje Z
# Equivale a inclinar la cabeza hacia un hombro (girar la cámara como un volante)
# En este caso, el eje Z (la dirección a la que miras) se mantiene fijo y la imagen gira puramente en el plano 2D

theta_z <- deg2rad(20)  # 20º

R_roll <- matrix(c(
    cos(theta_z), -sin(theta_z), 0,
    sin(theta_z),  cos(theta_z), 0,
    0,             0,            1
), nrow = 3, byrow = TRUE)


# Se pueden aplicar varios movimientos a la vez multiplicando las matrices
# En álgebra matricial, el orden de las multiplicaciones importa (la rotación no es conmutativa)
# Una convención común en cámaras es aplicar primero Roll, luego Pitch y finalmente Yaw
# lo que equivale a multiplicar en orden inverso de derecha a izquierda:
R_combinada <- R_yaw %*% R_pitch %*% R_roll

# Aplicar la transformación (Ej: tomada con lente de 28mm)
focal_length_mm <- 12
img_rotada <- rotate_camera(img, R_combinada, focal_length_mm, zoom = 0.5)
writeTIFF(img_rotada, "street12mm_rotate.tif")



#######################################
# 2. CORRECCIÓN DE TRAPEZOIDE A RELACIÓN DE ASPECTO REAL

# NOTA: con esta función de rotación, si en la escena aparece un objeto que en el mundo real es un rectángulo
# pero por efecto de las fugas (modelo de proyección pinhole) aparece como un trapezoide en la imagen origen
# y parametrizando una rotación logramos convertirlo en un rectángulo perfecto (lo que sería equivalente a
# apuntar la cámara perpendicularmente al plano que contiene dicho rectángulo), tendremos
# la relación de aspecto real del objeto, pero siempre y cuando la focal_length_mm fuera la exacta


# 1. Cargar una imagen de prueba
img <- readTIFF("building.tif")  # Fuji X-S10 con 11mm (16,5mm eq.)
# img=add_grid(img, n_gridx = 12, colour = c(1,1,0), linewidth = 4)
# writeTIFF(img, "building_grid.tif")

focal_length_mm <- 16.5


# Ejemplo: Las 4 esquinas de una fachada
# Orden: 1=Top-Left, 2=Bottom-Left, 3=Bottom-Right, 4=Top-Right
esquinas <- matrix(c(
    1516, 2391,   # p1: TL
    982, 3609,    # p2: BL
    6181, 3372,   # p3: BR
    5254, 1873    # p4: TR
), ncol=2, byrow=TRUE)

# 1. Obtener matriz rotación rectificadora
R_rect <- get_rectifying_rotation(H = nrow(img), W = ncol(img), fl_FF_mm = focal_length_mm, pts = esquinas)
img_frontal <- rotate_camera(img, R_rect, fl_FF_mm = focal_length_mm, zoom = 0.27)
writeTIFF(img_frontal, "building_frontal.tif")

# The resulting rectangle has dimensions: 2442 x 798 -> 3.06 aspect ratio, CORRRRRRRECT!!!
