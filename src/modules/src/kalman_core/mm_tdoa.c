/**
 * ,---------,       ____  _ __
 * |  ,-^-,  |      / __ )(_) /_______________ _____  ___
 * | (  O  ) |     / __  / / __/ ___/ ___/ __ `/_  / / _ \
 * | / ,--'  |    / /_/ / / /_/ /__/ /  / /_/ / / /_/  __/
 *    +------`   /_____/_/\__/\___/_/   \__,_/ /___/\___/
 *
 * Crazyflie control firmware
 *
 * Copyright (C) 2021 Bitcraze AB
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, in version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "mm_tdoa.h"
#include "outlierFilter.h"
#include "test_support.h"

// TODO krri What is this used for? Do we still need it?
TESTABLE_STATIC uint32_t tdoaCount = 0;

void kalmanCoreUpdateWithTDOA(kalmanCoreData_t* this, tdoaMeasurement_t *tdoa)
{
  if (tdoaCount >= 100)
  {
    /**
     * Measurement equation:
     * dR = dT + d1 - d0
     */

    float measurement = tdoa->distanceDiff;
    int iter_count = 0;
    float epsilon = 1e-10;
    float iter_error = 2*epsilon;
    kalmanCoreData_t *eta1 = this;
    kalmanCoreData_t *eta = eta1;
    while (iter_error > epsilon && iter_count < 25)
    {

      // predict based on current state
      float x = eta->S[KC_STATE_X];
      float y = eta->S[KC_STATE_Y];
      float z = eta->S[KC_STATE_Z];

      float x1 = tdoa->anchorPosition[1].x, y1 = tdoa->anchorPosition[1].y, z1 = tdoa->anchorPosition[1].z;
      float x0 = tdoa->anchorPosition[0].x, y0 = tdoa->anchorPosition[0].y, z0 = tdoa->anchorPosition[0].z;

      float dx1 = x - x1;
      float dy1 = y - y1;
      float dz1 = z - z1;

      float dy0 = y - y0;
      float dx0 = x - x0;
      float dz0 = z - z0;

      float d1 = sqrtf(powf(dx1, 2) + powf(dy1, 2) + powf(dz1, 2));
      float d0 = sqrtf(powf(dx0, 2) + powf(dy0, 2) + powf(dz0, 2));

      float predicted = d1 - d0;
      float error = measurement - predicted;

      float h[KC_STATE_DIM] = {0};
      arm_matrix_instance_f32 H = {1, KC_STATE_DIM, h};

      if ((d0 != 0.0f) && (d1 != 0.0f)) {
        h[KC_STATE_X] = (dx1 / d1 - dx0 / d0);
        h[KC_STATE_Y] = (dy1 / d1 - dy0 / d0);
        h[KC_STATE_Z] = (dz1 / d1 - dz0 / d0);

        vector_t jacobian = {
          .x = h[KC_STATE_X],
          .y = h[KC_STATE_Y],
          .z = h[KC_STATE_Z],
        };

        point_t estimatedPosition = {
          .x = eta->S[KC_STATE_X],
          .y = eta->S[KC_STATE_Y],
          .z = eta->S[KC_STATE_Z],
        };

        bool sampleIsGood = outlierFilterValidateTdoaSteps(tdoa, error, &jacobian, &estimatedPosition);
        if (sampleIsGood) {
          if(z > 1.0f){
            kalmanCoreScalarUpdate(eta1, &H, error, 0.8f*tdoa->stdDev);
          }
          else{
            kalmanCoreScalarUpdate(eta1, &H, error, tdoa->stdDev);
          }
        }
        
        /*eta1->S[KC_STATE_X] = this->S[KC_STATE_X] + eta->K_vector[KC_STATE_X] * (error - h[KC_STATE_X]*(this->S[KC_STATE_X]-eta->S[KC_STATE_X]));
        eta1->S[KC_STATE_Y] = this->S[KC_STATE_Y] + eta->K_vector[KC_STATE_Y] * (error - h[KC_STATE_Y]*(this->S[KC_STATE_Y]-eta->S[KC_STATE_Y]));
        eta1->S[KC_STATE_Z] = this->S[KC_STATE_Z] + eta->K_vector[KC_STATE_Z] * (error - h[KC_STATE_Z]*(this->S[KC_STATE_Z]-eta->S[KC_STATE_Z])); */

        float eta_norm = sqrtf(powf(eta->S[KC_STATE_X],2) + powf(eta->S[KC_STATE_Y],2) + powf(eta->S[KC_STATE_Z],2));
        float eta1_norm = sqrtf(powf(eta1->S[KC_STATE_X],2) + powf(eta1->S[KC_STATE_Y],2) + powf(eta1->S[KC_STATE_Z],2)); 
        iter_error = fabsf(eta1_norm - eta_norm)/ fabsf(eta_norm); 
      }
      eta = eta1;
      iter_count ++;
    }
    /*for (int i=0; i<KC_STATE_DIM; i++) {
      for (int j=i; j<KC_STATE_DIM; j++) {
        eta1->P[i][j] = eta->P[i][j];
        }} */
    this = eta1;
  }
  
  tdoaCount++;
}
