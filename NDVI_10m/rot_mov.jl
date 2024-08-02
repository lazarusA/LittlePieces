struct CameraRotationMove{T} <: FlyThroughPaths.PathChange{T}
    duration::T
    axis::Vec3{Float64}
    angle::Float64
    action
  end
  
  CameraRotationMove(duration::T, axis::Vec3, angle, action = nothing) where {T} = CameraRotationMove{T}(duration, axis, angle, action)
  
  function (move::CameraRotationMove)(view::ViewState, t)
  FlyThroughPaths.checkt(t, move)
  (; axis, angle, action) = move
  tf = t/FlyThroughPaths.duration(move)
  qrotation = Makie.qrotation(axis, angle * tf)
  viewing_direction = view.lookat - view.eyeposition
  rotation_matrix = Makie.rotationmatrix4(qrotation)[1:3, 1:3] # we don't need the fourth dimension
  new_view = ViewState(
    lookat = view.eyeposition + rotation_matrix * viewing_direction,
    eyeposition = view.eyeposition, upvector = view.upvector, fov = view.fov
  )
  new_view
  end