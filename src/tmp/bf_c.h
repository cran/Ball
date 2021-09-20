#ifndef BF_C_H
#define BF_C_H
 
    #ifdef __cplusplus
    extern "C" {
    #endif

    extern void *BF_C_new(double **distance_matrix, int fix_center_num, int num);
    extern void BF_C_delete(void *BF);
    extern void BF_C_train(void *BF);
    extern void BF_C_get_fitted(void *BF, double **predict_matrix);
    extern void BF_C_predict(void *BF, double **predict_matrix, double **new_distance, int pred_num);
    extern void BF_C_free_BF(void *BF);

    #ifdef __cplusplus
    }
    #endif
 
#endif
