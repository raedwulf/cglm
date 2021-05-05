/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Functions:
   CGLM_INLINE void  glm_frustum_lh_no(float left,    float right,
                                       float bottom,  float top,
                                       float nearVal, float farVal,
                                       mat4  dest)
   CGLM_INLINE void  glm_ortho_lh_no(float left,    float right,
                                       float bottom,  float top,
                                       float nearVal, float farVal,
                                       mat4  dest)
   CGLM_INLINE void  glm_perspective_lh_no(float fovy,
                                           float aspect,
                                           float nearVal,
                                           float farVal,
                                           mat4  dest)
   CGLM_INLINE void glm_ortho_aabb_lh_no(vec3 box[2], mat4 dest)
   CGLM_INLINE void glm_ortho_aabb_p_lh_no(vec3 box[2],
                                           float padding,
                                           mat4 dest)
   CGLM_INLINE void glm_ortho_aabb_pz_lh_no(vec3 box[2],
                                            float padding,
                                            mat4 dest)
   CGLM_INLINE void glm_ortho_default_lh_no(float aspect,
                                            mat4 dest)
   CGLM_INLINE void glm_ortho_default_s_lh_no(float aspect,
                                              float size,
                                              mat4 dest)
   CGLM_INLINE void glm_perspective_lh_no(float fovy,
                                          float aspect,
                                          float nearVal,
                                          float farVal,
                                          mat4  dest)
   CGLM_INLINE void glm_persp_move_far_lh_no(mat4 proj,
                                             float deltaFar)
   CGLM_INLINE void glm_persp_decomp_lh_no(mat4 proj,
                                           float * __restrict nearVal,
                                           float * __restrict farVal,
                                           float * __restrict top,
                                           float * __restrict bottom,
                                           float * __restrict left,
                                           float * __restrict right)
  CGLM_INLINE void glm_persp_decompv_lh_no(mat4 proj,
                                           float dest[6])
  CGLM_INLINE void glm_persp_decomp_x_lh_no(mat4 proj,
                                            float * __restrict left,
                                            float * __restrict right)
  CGLM_INLINE void glm_persp_decomp_y_lh_no(mat4 proj,
                                            float * __restrict top,
                                            float * __restrict bottom)
  CGLM_INLINE void glm_persp_decomp_z_lh_no(mat4 proj,
                                            float * __restrict nearVal,
                                            float * __restrict farVal)
  CGLM_INLINE void glm_persp_decomp_far_lh_no(mat4 proj, float * __restrict farVal)
  CGLM_INLINE void glm_persp_decomp_near_lh_no(mat4 proj, float * __restrict nearVal)
  CGLM_INLINE void glm_persp_sizes_lh_no(mat4 proj, float fovy, vec4 dest)
 */

#ifndef cglm_cam_lh_no_h
#define cglm_cam_lh_no_h

#include "common.h"
#include "plane.h"

/*!
 * @brief set up perspective peprojection matrix
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  left    viewport.left
 * @param[in]  right   viewport.right
 * @param[in]  bottom  viewport.bottom
 * @param[in]  top     viewport.top
 * @param[in]  nearVal near clipping plane
 * @param[in]  farVal  far clipping plane
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_frustum_lh_no(float left,    float right,
                  float bottom,  float top,
                  float nearVal, float farVal,
                  mat4  dest) {
  float rl, tb, fn, nv;

  glm_mat4_zero(dest);

  rl = 1.0f / (right  - left);
  tb = 1.0f / (top    - bottom);
  fn =-1.0f / (farVal - nearVal);
  nv = 2.0f * nearVal;

  dest[0][0] = nv * rl;
  dest[1][1] = nv * tb;
  dest[2][0] = (right  + left)    * rl;
  dest[2][1] = (top    + bottom)  * tb;
  dest[2][2] =-(farVal + nearVal) * fn;
  dest[2][3] = 1.0f;
  dest[3][2] = farVal * nv * fn;
}

/*!
 * @brief set up orthographic projection matrix
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  left    viewport.left
 * @param[in]  right   viewport.right
 * @param[in]  bottom  viewport.bottom
 * @param[in]  top     viewport.top
 * @param[in]  nearVal near clipping plane
 * @param[in]  farVal  far clipping plane
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_lh_no(float left,    float right,
                float bottom,  float top,
                float nearVal, float farVal,
                mat4  dest) {
  float rl, tb, fn;

  glm_mat4_zero(dest);

  rl = 1.0f / (right  - left);
  tb = 1.0f / (top    - bottom);
  fn =-1.0f / (farVal - nearVal);

  dest[0][0] = 2.0f * rl;
  dest[1][1] = 2.0f * tb;
  dest[2][2] =-2.0f * fn;
  dest[3][0] =-(right  + left)    * rl;
  dest[3][1] =-(top    + bottom)  * tb;
  dest[3][2] = (farVal + nearVal) * fn;
  dest[3][3] = 1.0f;
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box   AABB
 * @param[out] dest  result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_lh_no(vec3 box[2], mat4 dest) {
  glm_ortho_lh_no(box[0][0],  box[1][0],
                  box[0][1],  box[1][1],
                 -box[1][2], -box[0][2],
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_p_lh_no(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_lh_no(box[0][0] - padding,    box[1][0] + padding,
                  box[0][1] - padding,    box[1][1] + padding,
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding for near and far
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_pz_lh_no(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_lh_no(box[0][0],              box[1][0],
                  box[0][1],              box[1][1],
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up unit orthographic projection matrix
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  aspect aspect ration ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_ortho_default_lh_no(float aspect, mat4 dest) {
  if (aspect >= 1.0f) {
    glm_ortho_lh_no(-aspect, aspect, -1.0f, 1.0f, -100.0f, 100.0f, dest);
    return;
  }

  aspect = 1.0f / aspect;

  glm_ortho_lh_no(-1.0f, 1.0f, -aspect, aspect, -100.0f, 100.0f, dest);
}

/*!
 * @brief set up orthographic projection matrix with given CUBE size
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[in]  size   cube size
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_ortho_default_s_lh_no(float aspect, float size, mat4 dest) {
  if (aspect >= 1.0f) {
    glm_ortho_lh_no(-size * aspect,
                     size * aspect,
                    -size,
                     size,
                    -size - 100.0f,
                     size + 100.0f,
                     dest);
    return;
  }

  glm_ortho_lh_no(-size,
                   size,
                  -size / aspect,
                   size / aspect,
                  -size - 100.0f,
                   size + 100.0f,
                   dest);
}

/*!
 * @brief set up perspective projection matrix
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  fovy    field of view angle
 * @param[in]  aspect  aspect ratio ( width / height )
 * @param[in]  nearVal near clipping plane
 * @param[in]  farVal  far clipping planes
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_perspective_lh_no(float fovy,
                      float aspect,
                      float nearVal,
                      float farVal,
                      mat4  dest) {
  float f, fn;

  glm_mat4_zero(dest);

  f  = 1.0f / tanf(fovy * 0.5f);
  fn = 1.0f / (nearVal - farVal);

  dest[0][0] = f / aspect;
  dest[1][1] = f;
  dest[2][2] =-(nearVal + farVal) * fn;
  dest[2][3] = 1.0f;
  dest[3][2] = 2.0f * nearVal * farVal * fn;

}

/*!
 * @brief extend perspective projection matrix's far distance
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * this function does not guarantee far >= near, be aware of that!
 *
 * @param[in, out] proj      projection matrix to extend
 * @param[in]      deltaFar  distance from existing far (negative to shink)
 */
CGLM_INLINE
void
glm_persp_move_far_lh_no(mat4 proj, float deltaFar) {
  float fn, farVal, nearVal, p22, p32, h;

  p22        = -proj[2][2];
  p32        = proj[3][2];

  nearVal    = p32 / (p22 - 1.0f);
  farVal     = p32 / (p22 + 1.0f) + deltaFar;
  fn         = 1.0f / (nearVal - farVal);

  proj[2][2] = -(farVal + nearVal) * fn;
  proj[3][2] = 2.0f * nearVal * farVal * fn;
}

/*!
 * @brief decomposes frustum values of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 * @param[out] farVal  far
 * @param[out] top     top
 * @param[out] bottom  bottom
 * @param[out] left    left
 * @param[out] right   right
 */
CGLM_INLINE
void
glm_persp_decomp_lh_no(mat4 proj,
                       float * __restrict nearVal, float * __restrict farVal,
                       float * __restrict top,     float * __restrict bottom,
                       float * __restrict left,    float * __restrict right) {
  float m00, m11, m20, m21, m22, m32, n, f;
  float n_m11, n_m00;

  m00 = proj[0][0];
  m11 = proj[1][1];
  m20 = proj[2][0];
  m21 = proj[2][1];
  m22 =-proj[2][2];
  m32 = proj[3][2];

  n = m32 / (m22 - 1.0f);
  f = m32 / (m22 + 1.0f);

  n_m11 = n / m11;
  n_m00 = n / m00;

  *nearVal = n;
  *farVal  = f;
  *bottom  = n_m11 * (m21 - 1.0f);
  *top     = n_m11 * (m21 + 1.0f);
  *left    = n_m00 * (m20 - 1.0f);
  *right   = n_m00 * (m20 + 1.0f);
}

/*!
 * @brief decomposes frustum values of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *        this makes easy to get all values at once
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] dest   array
 */
CGLM_INLINE
void
glm_persp_decompv_lh_no(mat4 proj, float dest[6]) {
  glm_persp_decomp_lh_no(proj, &dest[0], &dest[1], &dest[2],
                               &dest[3], &dest[4], &dest[5]);
}

/*!
 * @brief decomposes left and right values of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *        x stands for x axis (left / right axis)
 *
 * @param[in]  proj  perspective projection matrix
 * @param[out] left  left
 * @param[out] right right
 */
CGLM_INLINE
void
glm_persp_decomp_x_lh_no(mat4 proj,
                         float * __restrict left,
                         float * __restrict right) {
  float nearVal, m20, m00;

  m00 = proj[0][0];
  m20 = proj[2][0];
  m22 =-proj[2][2];

  nearVal = proj[3][2] / (m22 - 1.0f);
  *left   = nearVal * (m20 - 1.0f) / m00;
  *right  = nearVal * (m20 + 1.0f) / m00;
}

/*!
 * @brief decomposes top and bottom values of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *        y stands for y axis (top / botom axis)
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] top    top
 * @param[out] bottom bottom
 */
CGLM_INLINE
void
glm_persp_decomp_y_lh_no(mat4 proj,
                         float * __restrict top,
                         float * __restrict bottom) {
  float nearVal, m21, m11, m22;

  m21 = proj[2][1];
  m11 = proj[1][1];
  m22 =-proj[2][2];

  nearVal = proj[3][2] / (m22 - 1.0f);
  *bottom = nearVal * (m21 - 1.0f) / m11;
  *top    = nearVal * (m21 + 1.0f) / m11;
}

/*!
 * @brief decomposes near and far values of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *        z stands for z axis (near / far axis)
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 * @param[out] farVal  far
 */
CGLM_INLINE
void
glm_persp_decomp_z_lh_no(mat4 proj,
                      float * __restrict nearVal,
                      float * __restrict farVal) {
  float m32, m22;

  m32 = proj[3][2];
  m22 =-proj[2][2];

  *nearVal = m32 / (m22 - 1.0f);
  *farVal  = m32 / (m22 + 1.0f);
}

/*!
 * @brief decomposes far value of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] farVal far
 */
CGLM_INLINE
void
glm_persp_decomp_far_lh_no(mat4 proj, float * __restrict farVal) {
  *farVal = proj[3][2] / (-proj[2][2] + 1.0f);
}

/*!
 * @brief decomposes near value of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 */
CGLM_INLINE
void
glm_persp_decomp_near_lh_no(mat4 proj, float * __restrict nearVal) {
  *nearVal = proj[3][2] / (-proj[2][2] - 1.0f);
}

/*!
 * @brief returns sizes of near and far planes of perspective projection
 *        with a left-hand coordinate system and a
 *        clip-space of [-1, 1].
 *
 * @param[in]  proj perspective projection matrix
 * @param[in]  fovy fovy (see brief)
 * @param[out] dest sizes order: [Wnear, Hnear, Wfar, Hfar]
 */
CGLM_INLINE
void
glm_persp_sizes_lh_no(mat4 proj, float fovy, vec4 dest) {
  float t, a, nearVal, farVal;

  t = 2.0f * tanf(fovy * 0.5f);
  a = glm_persp_aspect(proj);

  glm_persp_decomp_z_lh_no(proj, &nearVal, &farVal);

  dest[1]  = t * nearVal;
  dest[3]  = t * farVal;
  dest[0]  = a * dest[1];
  dest[2]  = a * dest[3];
}

#endif /*cglm_cam_lh_no_h*/
