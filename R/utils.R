
parFun <- function(m1,m2, 
                   bx,
                   bm1x,
                   bm2x,
                   bm1m2x) {
  bx + m1 * bm1x  + m2 * bm2x  + m1 * m2 * bm1m2x 
}

seFun <- function(m1,m2,
                  Vx,
                  Vm1x,
                  Vm2x,
                  Vm1m2x,
                  covX_m1x,
                  covX_m2x,
                  covX_m1m2x,
                  covm1x_m2x,
                  covm1x_m1m2x,
                  covm2x_m1m2x) {
  sqrt(Vx +  
         m1^2 * Vm1x +
         m2^2 * Vm2x +
         (m1 * m2)^2 * Vm1m2x +
         2 * m1 * covX_m1x +
         2 * m2 * covX_m2x +
         2 * m1 * m2 * covX_m1m2x +
         2 * m1 * m2 * covm1x_m2x +
         2 * m1^2 * m2 * covm1x_m1m2x +
         2 * m1 * m2^2 * covm2x_m1m2x)
}
