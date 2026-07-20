#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rotate_camera_rcpp(NumericVector img, NumericMatrix R, double fl_FF_mm, double zoom = 1.0, int shift_x = 0, int shift_y = 0) {
    // 1. Extraer dimensiones del array
    IntegerVector dims = img.attr("dim");
    if (dims.size() != 3) {
        stop("La imagen debe ser un array de 3 dimensiones (H x W x 3).");
    }
    int H = dims[0];
    int W = dims[1];
    int C = dims[2];

    // 2. Crear array de salida inicializado con ceros
    NumericVector out(img.size());
    out.attr("dim") = dims;
    
    // 3. Conversión robusta de distancia focal (diagonal Full Frame = 43.2666 mm)
    double diag_mm = std::sqrt(36.0 * 36.0 + 24.0 * 24.0); 
    double diag_pixel = std::sqrt((double)H * H + (double)W * W);
    
    // Distancia focal de la imagen original (para de-proyectar al 3D)
    double f_pixel_src = fl_FF_mm * (diag_pixel / diag_mm);
    
    // Distancia focal de la cámara destino modificada por el zoom
    // zoom > 1 reduce el campo de visión (acerca)
    // zoom < 1 amplía el campo de visión (aleja)
    double f_pixel_dst = f_pixel_src * zoom;

    // Centro óptico
    double cx = (W - 1.0) / 2.0;
    double cy = (H - 1.0) / 2.0;

    // 4. Matriz de rotación inversa (transpuesta)
    double R_inv[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            R_inv[i][j] = R(j, i); 
        }
    }

    int plane_size = H * W; 

    // 5. Bucle principal de mapeo inverso
    for (int x = 0; x < W; ++x) {
        for (int y = 0; y < H; ++y) {
            
            // Paso A: De-proyección (DESTINO). 
            // Aplicamos el desplazamiento shift_x y shift_y al mapeo inverso
            double ray_x = (x - shift_x - cx) / f_pixel_dst;
            double ray_y = (y - shift_y - cy) / f_pixel_dst;
            double ray_z = 1.0;

            // Paso B: Rotar el rayo hacia el origen
            double src_ray_x = R_inv[0][0]*ray_x + R_inv[0][1]*ray_y + R_inv[0][2]*ray_z;
            double src_ray_y = R_inv[1][0]*ray_x + R_inv[1][1]*ray_y + R_inv[1][2]*ray_z;
            double src_ray_z = R_inv[2][0]*ray_x + R_inv[2][1]*ray_y + R_inv[2][2]*ray_z;

            // Descartar rayos que apuntan hacia atrás
            if (src_ray_z <= 0.0) continue;

            // Paso C: Proyección (ORIGEN). Usamos f_pixel_src inalterada
            double src_x = f_pixel_src * (src_ray_x / src_ray_z) + cx;
            double src_y = f_pixel_src * (src_ray_y / src_ray_z) + cy;

            // Paso D: Interpolación bilineal
            if (src_x >= 0.0 && src_x < (W - 1) && src_y >= 0.0 && src_y < (H - 1)) {
                
                int x0 = std::floor(src_x);
                int y0 = std::floor(src_y);
                int x1 = x0 + 1;
                int y1 = y0 + 1;

                double dx = src_x - x0;
                double dy = src_y - y0;

                double w00 = (1.0 - dx) * (1.0 - dy);
                double w10 = dx * (1.0 - dy);
                double w01 = (1.0 - dx) * dy;
                double w11 = dx * dy;

                for (int c = 0; c < C; ++c) {
                    int c_offset = c * plane_size;
                    
                    double val = w00 * img[y0 + x0 * H + c_offset] +
                                 w10 * img[y0 + x1 * H + c_offset] +
                                 w01 * img[y1 + x0 * H + c_offset] +
                                 w11 * img[y1 + x1 * H + c_offset];

                    out[y + x * H + c_offset] = val;
                }
            }
        }
    }

    return out;
}
