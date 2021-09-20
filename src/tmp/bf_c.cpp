#include "bf_c.h"
#include "bf.h"
 
extern void *BF_C_new(double** distance_matrix, int fix_center_num, int num) {
    return new BF(distance_matrix, fix_center_num, num);
}
 
extern void BF_C_delete(void *bf_object) {
    BF *bf = (BF *) bf_object;
    delete bf;
}
 
extern void BF_C_train(void *bf_object) {
    BF *bf = (BF *) bf_object;
    return bf->train();
}

extern void BF_C_get_fitted(void *bf_object, double **predict_matrix) {
    BF *bf = (BF *) bf_object;
    return bf->get_fitted(predict_matrix);
}

extern void BF_C_predict(void *bf_object, double **predict_matrix, double **new_distance, int pred_num) {
    BF *bf = (BF *) bf_object;
    return bf->predict(predict_matrix, new_distance, pred_num);
}

extern void BF_C_free_BF(void *bf_object) {
    BF *bf = (BF *) bf_object;
    return bf->free_BF();
}
