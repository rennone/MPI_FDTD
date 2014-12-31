#ifndef _EVALATE_H__
#define _EVALATE_H__

enum EvalKinds
{
  EVAL_RED,
  EVAL_GREEN,
  EVAL_BLUE,
  EVAL_NUM
};

// TODO : EvalKindsを引数にとるようにするまで使ってはいけない
//各波長の角反射率
// reflec[lambda][degree] : 波長λのdegree°方向の反射率
extern double evaluate_evaluate(double **reflec, int stLambda, int enLambda);

extern double evaluate_eval(double h, double s, double v, enum EvalKinds k );

#endif
