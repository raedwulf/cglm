/*
 * Copyright (c), Recep Aslantas.
 *
 * MIT License (MIT), http://opensource.org/licenses/MIT
 * Full license can be found in the LICENSE file
 */

/*
 Functions:
   CGLM_INLINE void  glm_frustum(float left,    float right,
                                 float bottom,  float top,
                                 float nearVal, float farVal,
                                 mat4  dest)
   CGLM_INLINE void  glm_ortho(float left,    float right,
                               float bottom,  float top,
                               float nearVal, float farVal,
                               mat4  dest)
   CGLM_INLINE void  glm_ortho_aabb(vec3 box[2], mat4 dest)
   CGLM_INLINE void  glm_ortho_aabb_p(vec3 box[2],  float padding, mat4 dest)
   CGLM_INLINE void  glm_ortho_aabb_pz(vec3 box[2], float padding, mat4 dest)
   CGLM_INLINE void  glm_ortho_default(float aspect, mat4  dest)
   CGLM_INLINE void  glm_ortho_default_s(float aspect, float size, mat4 dest)
   CGLM_INLINE void  glm_perspective(float fovy,
                                     float aspect,
                                     float nearVal,
                                     float farVal,
                                     mat4  dest)
   CGLM_INLINE void  glm_perspective_default(float aspect, mat4 dest)
   CGLM_INLINE void  glm_perspective_resize(float aspect, mat4 proj)
   CGLM_INLINE void  glm_lookat(vec3 eye, vec3 center, vec3 up, mat4 dest)
   CGLM_INLINE void  glm_look(vec3 eye, vec3 dir, vec3 up, mat4 dest)
   CGLM_INLINE void  glm_look_anyup(vec3 eye, vec3 dir, mat4 dest)
   CGLM_INLINE void  glm_persp_decomp(mat4   proj,
                                      float *nearVal, float *farVal,
                                      float *top,     float *bottom,
                                      float *left,    float *right)
   CGLM_INLINE void  glm_persp_decompv(mat4 proj, float dest[6])
   CGLM_INLINE void  glm_persp_decomp_x(mat4 proj, float *left, float *right)
   CGLM_INLINE void  glm_persp_decomp_y(mat4 proj, float *top,  float *bottom)
   CGLM_INLINE void  glm_persp_decomp_z(mat4 proj, float *nearv, float *farv)
   CGLM_INLINE void  glm_persp_decomp_far(mat4 proj, float *farVal)
   CGLM_INLINE void  glm_persp_decomp_near(mat4 proj, float *nearVal)
   CGLM_INLINE float glm_persp_fovy(mat4 proj)
   CGLM_INLINE float glm_persp_aspect(mat4 proj)
   CGLM_INLINE void  glm_persp_sizes(mat4 proj, float fovy, vec4 dest)
 */

#ifndef cglm_vcam_h
#define cglm_vcam_h

#include "common.h"
#include "plane.h"

/*!
 * @brief set up perspective peprojection matrix
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
glm_frustum(float left,    float right,
            float bottom,  float top,
            float nearVal, float farVal,
            mat4  dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_frustum_lh_zo(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_frustum_lh_no(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_frustum_rh_zo(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_frustum_rh_no(left, right, bottom, top, nearVal, farVal, dest);
#endif
}

/*!
 * @brief set up orthographic projection matrix
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
glm_ortho(float left,    float right,
          float bottom,  float top,
          float nearVal, float farVal,
          mat4  dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_ortho_lh_zo(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_ortho_lh_no(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_ortho_rh_zo(left, right, bottom, top, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_ortho_rh_no(left, right, bottom, top, nearVal, farVal, dest);
#endif
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & ZO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box   AABB
 * @param[out] dest  result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_rh_zo(vec3 box[2], mat4 dest) {
  glm_ortho_rh_zo(box[0][0],  box[1][0],
                  box[0][1],  box[1][1],
                 -box[1][2], -box[0][2],
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & NO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box   AABB
 * @param[out] dest  result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_rh_no(vec3 box[2], mat4 dest) {
  glm_ortho_rh_no(box[0][0],  box[1][0],
                  box[0][1],  box[1][1],
                 -box[1][2], -box[0][2],
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box   AABB
 * @param[out] dest  result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb(vec3 box[2], mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_ortho_lh_zo(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_ortho_lh_no(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_ortho_rh_zo(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_ortho_rh_no(box, dest);
#endif
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & ZO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_p_rh_zo(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_rh_zo(box[0][0] - padding,    box[1][0] + padding,
                  box[0][1] - padding,    box[1][1] + padding,
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & NO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_p_rh_no(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_rh_no(box[0][0] - padding,    box[1][0] + padding,
                  box[0][1] - padding,    box[1][1] + padding,
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_p(vec3 box[2], float padding, mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_ortho_aabb_p_lh_zo(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_ortho_aabb_p_lh_no(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_ortho_aabb_p_rh_zo(box, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_ortho_aabb_p_rh_no(box, dest);
#endif
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & ZO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding for near and far
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_pz_rh_zo(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_rh_zo(box[0][0],              box[1][0],
                  box[0][1],              box[1][1],
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box (RH & NO)
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding for near and far
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_pz_rh_no(vec3 box[2], float padding, mat4 dest) {
  glm_ortho_rh_no(box[0][0],              box[1][0],
                  box[0][1],              box[1][1],
                -(box[1][2] + padding), -(box[0][2] - padding),
                  dest);
}

/*!
 * @brief set up orthographic projection matrix using bounding box
 *
 * bounding box (AABB) must be in view space
 *
 * @param[in]  box     AABB
 * @param[in]  padding padding for near and far
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_ortho_aabb_pz(vec3 box[2], float padding, mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_ortho_aabb_pz_lh_zo(box, padding, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_ortho_aabb_pz_lh_no(box, padding, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_ortho_aabb_pz_rh_zo(box, padding, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_ortho_aabb_pz_rh_no(box, padding, dest);
#endif
}

/*!
 * @brief set up unit orthographic projection matrix
 *
 * @param[in]  aspect aspect ration ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_ortho_default(float aspect, mat4 dest) {
  if (aspect >= 1.0f) {
    glm_ortho(-aspect, aspect, -1.0f, 1.0f, -100.0f, 100.0f, dest);
    return;
  }

  aspect = 1.0f / aspect;

  glm_ortho(-1.0f, 1.0f, -aspect, aspect, -100.0f, 100.0f, dest);
}

/*!
 * @brief set up orthographic projection matrix with given CUBE size
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[in]  size   cube size
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_ortho_default_s(float aspect, float size, mat4 dest) {
  if (aspect >= 1.0f) {
    glm_ortho(-size * aspect,
               size * aspect,
              -size,
               size,
              -size - 100.0f,
               size + 100.0f,
               dest);
    return;
  }

  glm_ortho(-size,
             size,
            -size / aspect,
             size / aspect,
            -size - 100.0f,
             size + 100.0f,
             dest);
}

/*!
 * @brief set up perspective projection matrix
 *
 * @param[in]  fovy    field of view angle
 * @param[in]  aspect  aspect ratio ( width / height )
 * @param[in]  nearVal near clipping plane
 * @param[in]  farVal  far clipping planes
 * @param[out] dest    result matrix
 */
CGLM_INLINE
void
glm_perspective(float fovy,
                float aspect,
                float nearVal,
                float farVal,
                mat4  dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_perspective_lh_zo(fovy, aspect, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_perspective_lh_no(fovy, aspect, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_perspective_rh_zo(fovy, aspect, nearVal, farVal, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_perspective_rh_no(fovy, aspect, nearVal, farVal, dest);
#endif
}

/*!
 * @brief extend perspective projection matrix's far distance
 *
 * this function does not guarantee far >= near, be aware of that!
 *
 * @param[in, out] proj      projection matrix to extend
 * @param[in]      deltaFar  distance from existing far (negative to shink)
 */
CGLM_INLINE
void
glm_persp_move_far(mat4 proj, float deltaFar) {
  float fn, farVal, nearVal, p22, p32, h;

  h          =-proj[2][3];
  p22        = proj[2][2] * h;
  p32        = proj[3][2];

  nearVal    = p32 / (p22 - 1.0f);
  farVal     = p32 / (p22 + 1.0f) + deltaFar;
  fn         = 1.0f / (nearVal - farVal);

  proj[2][2] = (nearVal + farVal) * fn;
  proj[3][2] = 2.0f * nearVal * farVal * fn;
}

/*!
 * @brief set up perspective projection matrix with default near/far
 *        and angle values (LH & NO)
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_perspective_default_lh_no(float aspect, mat4 dest) {
  glm_perspective_lh_no(GLM_PI_4f, aspect, 0.01f, 100.0f, dest);
}

/*!
 * @brief set up perspective projection matrix with default near/far
 *        and angle values (RH & ZO)
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_perspective_default_rh_zo(float aspect, mat4 dest) {
  glm_perspective_rh_zo(GLM_PI_4f, aspect, 0.01f, 100.0f, dest);
}

/*!
 * @brief set up perspective projection matrix with default near/far
 *        and angle values (RH & NO)
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_perspective_default_rh_no(float aspect, mat4 dest) {
  glm_perspective_rh_no(GLM_PI_4f, aspect, 0.01f, 100.0f, dest);
}

/*!
 * @brief set up perspective projection matrix with default near/far
 *        and angle values
 *
 * @param[in]  aspect aspect ratio ( width / height )
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_perspective_default(float aspect, mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_ZO
  glm_perspective_default_lh_zo(aspect, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH_NO
  glm_perspective_default_lh_no(aspect, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_ZO
  glm_perspective_default_rh_zo(aspect, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH_NO
  glm_perspective_default_rh_no(aspect, dest);
#endif

/*!
 * @brief resize perspective matrix by aspect ratio ( width / height )
 *        this makes very easy to resize proj matrix when window /viewport
 *        reized
 *
 * @param[in]      aspect aspect ratio ( width / height )
 * @param[in, out] proj   perspective projection matrix
 */
CGLM_INLINE
void
glm_perspective_resize(float aspect, mat4 proj) {
  if (proj[0][0] == 0.0f)
    return;

  proj[0][0] = proj[1][1] / aspect;
}

/*!
 * @brief set up view matrix
 *
 * NOTE: The UP vector must not be parallel to the line of sight from
 *       the eye point to the reference point
 *
 * @param[in]  eye    eye vector
 * @param[in]  center center vector
 * @param[in]  up     up vector
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_lookat(vec3 eye, vec3 center, vec3 up, mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH
  glm_lookat_lh(eye, center, up, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH
  glm_lookat_rh(eye, center, up, dest);
#endif
}

/*!
 * @brief set up view matrix
 *
 * convenient wrapper for lookat: if you only have direction not target self
 * then this might be useful. Because you need to get target from direction.
 *
 * NOTE: The UP vector must not be parallel to the line of sight from
 *       the eye point to the reference point
 *
 * @param[in]  eye    eye vector
 * @param[in]  dir    direction vector
 * @param[in]  up     up vector
 * @param[out] dest   result matrix
 */
CGLM_INLINE
void
glm_look(vec3 eye, vec3 dir, vec3 up, mat4 dest) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_LH
  glm_look_lh(eye, dir, up, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_RH
  glm_look_rh(eye, dir, up, dest);
#endif
}

/*!
 * @brief decomposes frustum values of perspective projection. (ZO)
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
glm_persp_decomp_zo(mat4 proj,
                    float * __restrict nearVal, float * __restrict farVal,
                    float * __restrict top,     float * __restrict bottom,
                    float * __restrict left,    float * __restrict right) {
  float m00, m11, m20, m21, m22, m32, n, f;
  float n_m11, n_m00;

  m00 = proj[0][0];
  m11 = proj[1][1];
  m20 = proj[2][0];
  m21 = proj[2][1];
  m22 = proj[2][2] * -proj[2][3];
  m32 = proj[3][2];

  n = m32 / m22;
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
 * @brief decomposes frustum values of perspective projection. (NO)
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
glm_persp_decomp_no(mat4 proj,
                    float * __restrict nearVal, float * __restrict farVal,
                    float * __restrict top,     float * __restrict bottom,
                    float * __restrict left,    float * __restrict right) {
  float m00, m11, m20, m21, m22, m32, n, f;
  float n_m11, n_m00;

  m00 = proj[0][0];
  m11 = proj[1][1];
  m20 = proj[2][0];
  m21 = proj[2][1];
  m22 = proj[2][2] * -proj[2][3];
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
 * @brief decomposes frustum values of perspective projection.
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
glm_persp_decomp(mat4 proj,
                 float * __restrict nearVal, float * __restrict farVal,
                 float * __restrict top,     float * __restrict bottom,
                 float * __restrict left,    float * __restrict right) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decomp_zo(proj, nearVal, farVal, top, bottom, left, right);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decomp_no(proj, nearVal, farVal, top, bottom, left, right);
#endif
}

/*!
 * @brief decomposes frustum values of perspective projection. (ZO)
 *        this makes easy to get all values at once
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] dest   array
 */
CGLM_INLINE
void
glm_persp_decompv_zo(mat4 proj, float dest[6]) {
  glm_persp_decomp_zo(proj, &dest[0], &dest[1], &dest[2],
                            &dest[3], &dest[4], &dest[5]);
}

/*!
 * @brief decomposes frustum values of perspective projection. (NO)
 *        this makes easy to get all values at once
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] dest   array
 */
CGLM_INLINE
void
glm_persp_decompv_no(mat4 proj, float dest[6]) {
  glm_persp_decomp_no(proj, &dest[0], &dest[1], &dest[2],
                            &dest[3], &dest[4], &dest[5]);
}

/*!
 * @brief decomposes frustum values of perspective projection.
 *        this makes easy to get all values at once
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] dest   array
 */
CGLM_INLINE
void
glm_persp_decompv(mat4 proj, float dest[6]) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decompv_zo(proj, dest);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decompv_no(proj, dest);
#endif
}

/*!
 * @brief decomposes left and right values of perspective projection (ZO).
 *        x stands for x axis (left / right axis)
 *
 * @param[in]  proj  perspective projection matrix
 * @param[out] left  left
 * @param[out] right right
 */
CGLM_INLINE
void
glm_persp_decomp_x_zo(mat4 proj,
                      float * __restrict left,
                      float * __restrict right) {
  float nearVal, m20, m00;

  m00 = proj[0][0];
  m20 = proj[2][0];

  nearVal = proj[3][2] / (proj[3][3]);
  *left   = nearVal * (m20 - 1.0f) / m00;
  *right  = nearVal * (m20 + 1.0f) / m00;
}

/*!
 * @brief decomposes left and right values of perspective projection (NO).
 *        x stands for x axis (left / right axis)
 *
 * @param[in]  proj  perspective projection matrix
 * @param[out] left  left
 * @param[out] right right
 */
CGLM_INLINE
void
glm_persp_decomp_x_no(mat4 proj,
                      float * __restrict left,
                      float * __restrict right) {
  float nearVal, m20, m00;

  m00 = proj[0][0];
  m20 = proj[2][0];

  nearVal = proj[3][2] / (proj[3][3] - 1.0f);
  *left   = nearVal * (m20 - 1.0f) / m00;
  *right  = nearVal * (m20 + 1.0f) / m00;
}

/*!
 * @brief decomposes left and right values of perspective projection.
 *        x stands for x axis (left / right axis)
 *
 * @param[in]  proj  perspective projection matrix
 * @param[out] left  left
 * @param[out] right right
 */
CGLM_INLINE
void
glm_persp_decomp_x(mat4 proj,
                   float * __restrict left,
                   float * __restrict right) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decomp_x_zo(proj, left, right);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decomp_x_no(proj, left, right);
#endif
}

/*!
 * @brief decomposes top and bottom values of perspective projection (ZO).
 *        y stands for y axis (top / botom axis)
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] top    top
 * @param[out] bottom bottom
 */
CGLM_INLINE
void
glm_persp_decomp_y_zo(mat4 proj,
                      float * __restrict top,
                      float * __restrict bottom) {
  float nearVal, m21, m11;

  m21 = proj[2][1];
  m11 = proj[1][1];

  nearVal = proj[3][2] / (proj[3][3]);
  *bottom = nearVal * (m21 - 1) / m11;
  *top    = nearVal * (m21 + 1) / m11;
}

/*!
 * @brief decomposes top and bottom values of perspective projection (NO).
 *        y stands for y axis (top / botom axis)
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] top    top
 * @param[out] bottom bottom
 */
CGLM_INLINE
void
glm_persp_decomp_y_no(mat4 proj,
                      float * __restrict top,
                      float * __restrict bottom) {
  float nearVal, m21, m11;

  m21 = proj[2][1];
  m11 = proj[1][1];

  nearVal = proj[3][2] / (proj[3][3] - 1.0f);
  *bottom = nearVal * (m21 - 1) / m11;
  *top    = nearVal * (m21 + 1) / m11;
}

/*!
 * @brief decomposes top and bottom values of perspective projection.
 *        y stands for y axis (top / botom axis)
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] top    top
 * @param[out] bottom bottom
 */
CGLM_INLINE
void
glm_persp_decomp_y(mat4 proj,
                   float * __restrict top,
                   float * __restrict bottom) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decomp_y_zo(proj, top, bottom);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decomp_y_no(proj, top, bottom);
#endif
}

/*!
 * @brief decomposes near and far values of perspective projection (ZO).
 *        z stands for z axis (near / far axis)
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 * @param[out] farVal  far
 */
CGLM_INLINE
void
glm_persp_decomp_z_zo(mat4 proj,
                      float * __restrict nearVal,
                      float * __restrict farVal) {
  float m32, m22;

  m32 = proj[3][2];
  m22 = proj[2][2] * -proj[2][3];

  *nearVal = m32 / m22;
  *farVal  = m32 / (m22 + 1.0f);
}

/*!
 * @brief decomposes near and far values of perspective projection (NO).
 *        z stands for z axis (near / far axis)
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 * @param[out] farVal  far
 */
CGLM_INLINE
void
glm_persp_decomp_z_no(mat4 proj,
                      float * __restrict nearVal,
                      float * __restrict farVal) {
  float m32, m22;

  m32 = proj[3][2];
  m22 = proj[2][2] * -proj[2][3];

  *nearVal = m32 / (m22 - 1.0f);
  *farVal  = m32 / (m22 + 1.0f);
}

/*!
 * @brief decomposes near and far values of perspective projection.
 *        z stands for z axis (near / far axis)
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 * @param[out] farVal  far
 */
CGLM_INLINE
void
glm_persp_decomp_z(mat4 proj,
                   float * __restrict nearVal,
                   float * __restrict farVal) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decomp_z_zo(proj, nearVal, farVal);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decomp_z_no(proj, nearVal, farVal);
#endif
}

/*!
 * @brief decomposes far value of perspective projection.
 *
 * @param[in]  proj   perspective projection matrix
 * @param[out] farVal far
 */
CGLM_INLINE
void
glm_persp_decomp_far(mat4 proj, float * __restrict farVal) {
  *farVal = proj[3][2] / (proj[2][2] * -proj[2][3] + 1.0f);
}

/*!
 * @brief decomposes near value of perspective projection (ZO).
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 */
CGLM_INLINE
void
glm_persp_decomp_near_zo(mat4 proj, float * __restrict nearVal) {
  *nearVal = proj[3][2] / (proj[2][2] * -proj[2][3]);
}

/*!
 * @brief decomposes near value of perspective projection (NO).
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 */
CGLM_INLINE
void
glm_persp_decomp_near_no(mat4 proj, float * __restrict nearVal) {
  *nearVal = proj[3][2] / (proj[2][2] * -proj[2][3] - 1.0f);
}

/*!
 * @brief decomposes near value of perspective projection.
 *
 * @param[in]  proj    perspective projection matrix
 * @param[out] nearVal near
 */
CGLM_INLINE
void
glm_persp_decomp_near(mat4 proj, float * __restrict nearVal) {
#if CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_ZO
  glm_persp_decomp_near_zo(proj, nearVal);
#elif CGLM_CONFIG_CLIP_CONTROL & CGLM_CLIP_CONTROL_NO
  glm_persp_decomp_near_no(proj, nearVal);
#endif
}

/*!
 * @brief returns field of view angle along the Y-axis (in radians)
 *
 * if you need to degrees, use glm_deg to convert it or use this:
 * fovy_deg = glm_deg(glm_persp_fovy(projMatrix))
 *
 * @param[in] proj perspective projection matrix
 */
CGLM_INLINE
float
glm_persp_fovy(mat4 proj) {
  return 2.0f * atanf(1.0f / proj[1][1]);
}

/*!
 * @brief returns aspect ratio of perspective projection
 *
 * @param[in] proj perspective projection matrix
 */
CGLM_INLINE
float
glm_persp_aspect(mat4 proj) {
  return proj[1][1] / proj[0][0];
}

/*!
 * @brief returns sizes of near and far planes of perspective projection
 *
 * @param[in]  proj perspective projection matrix
 * @param[in]  fovy fovy (see brief)
 * @param[out] dest sizes order: [Wnear, Hnear, Wfar, Hfar]
 */
CGLM_INLINE
void
glm_persp_sizes(mat4 proj, float fovy, vec4 dest) {
  float t, a, nearVal, farVal;

  t = 2.0f * tanf(fovy * 0.5f);
  a = glm_persp_aspect(proj);

  glm_persp_decomp_z(proj, &nearVal, &farVal);

  dest[1]  = t * nearVal;
  dest[3]  = t * farVal;
  dest[0]  = a * dest[1];
  dest[2]  = a * dest[3];
}

#endif /* cglm_vcam_h */
