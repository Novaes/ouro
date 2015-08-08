/*******************************************************************************
!  Copyright(C) 2001-2015 Intel Corporation. All Rights Reserved.
!
!  The source code, information  and  material ("Material") contained herein is
!  owned  by Intel Corporation or its suppliers or licensors, and title to such
!  Material remains  with Intel Corporation  or its suppliers or licensors. The
!  Material  contains proprietary information  of  Intel or  its  suppliers and
!  licensors. The  Material is protected by worldwide copyright laws and treaty
!  provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!  modified, published, uploaded, posted, transmitted, distributed or disclosed
!  in any way  without Intel's  prior  express written  permission. No  license
!  under  any patent, copyright  or  other intellectual property rights  in the
!  Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!  implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!  intellectual  property  rights must  be express  and  approved  by  Intel in
!  writing.
!
!  *Third Party trademarks are the property of their respective owners.
!
!  Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!  this  notice or  any other notice embedded  in Materials by Intel or Intel's
!  suppliers or licensors in any way.
!
!*******************************************************************************
!  Content:
!  Example of 2-dimension linear convolution operation on double precision data.
!*******************************************************************************/

#include <math.h>
#include <stdio.h>

#include "mkl_vsl.h"

int main()
{
    VSLConvTaskPtr task;
    double x[3*2]={1,1,1, 1,1,1};
    double y[4*3]={1,1,1,1, 1,1,1,1, 1,1,1,1};
    double z[6*4]={0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0};
    double e[6*4]={1,2,3,3,2,1, 2,4,6,6,4,2, 2,4,6,6,4,2, 1,2,3,3,2,1};
    MKL_INT xshape[2]={3,2}, yshape[2]={4,3}, zshape[2]={6,4};
    MKL_INT rank=2;
    int status,i,j;

    int mode = VSL_CONV_MODE_AUTO;

    /*
    *  Create task descriptor (create descriptor of problem)
    */
    status = vsldConvNewTask(&task,mode,rank,xshape,yshape,zshape);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: creation of job failed, exit with %d\n", status);
        return 1;
    }

    /*
    *  Execute task (Calculate 2 dimension convolution of two arrays)
    */
    status = vsldConvExec(task,x,NULL,y,NULL,z,NULL);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: job status bad, exit with %d\n", status);
        return 1;
    }

    /*
    *  Delete task object (delete descriptor of problem)
    */
    status = vslConvDeleteTask(&task);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: failed to delete task object, exit with %d\n", status);
        return 1;
    }

    /*
    * Check resulst for correctness:
    */
    for (j=0; j<zshape[1]; j++)
        for (i=0; i<zshape[0]; i++) {
            double zij = z[i + zshape[0]*j];
            double eij = e[i + zshape[0]*j];
            if (fabs(zij-eij) > fabs(eij)*1e-10) {
                printf("ERROR: wrong results:\n");
                printf("    z[%2d,%2d]: %lg\n",i,j,zij);
                printf("    expected: %lg\n",eij);
                printf("EXAMPLE FAILED\n");
                return 1;
            }
        }

    printf("EXAMPLE PASSED\n");
    return 0;
}