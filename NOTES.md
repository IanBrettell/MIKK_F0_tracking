# Notes

## 20220211

Noticed errors with tracking:

* `/nfs/research/birney/projects/indigene/MIKK_F0_tracking/split/novel_object/20191112_1236_18-2_L_A_q1.avi`
* `/nfs/research/birney/projects/indigene/MIKK_F0_tracking/split/novel_object/20191114_1305_50-2_R_A_q4.avi`

Caused by fish remaining still for most of the video, which created "ghost" segments. 

Fixed by setting bgsub = 'False' and adjusting intensity and area thresholds.

Code is configured in a way that only allows you to change parameters for all four quadrants. 

## 20220203

`/nfs/research/birney/users/ian/MIKK_F0_tracking/split/novel_object/20191119_1100_139-4_L_B_q4.mp4` takes days to track due to repeated CNN retraining that fails to reach the limit ratio of images to be used during pretraining: (e.g. `0.8090 (if higher than 0.95 we stop)`

When viewing the video with these parameters (Intensity `[0,200]`, Area `[150, 60000]`), it fails to pick up on one of the fishes during part of the video, and creates a "shadow" segment for one or both fishes.
This may be caused by the fact that neither fish moves for half of the video.
Try tracking with Intensity `[0,220]` which prevents the loss of segment for the second (MIKK) fish, but creates more "shadows".

View tracked video to confirm result.
