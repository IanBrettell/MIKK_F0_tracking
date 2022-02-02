# Notes

## 20220203

`/nfs/research/birney/users/ian/MIKK_F0_tracking/split/novel_object/20191119_1100_139-4_L_B_q4.mp4` takes days to track due to repeated CNN retraining that fails to reach the limit ratio of images to be used during pretraining: (e.g. `0.8090 (if higher than 0.95 we stop)`

When viewing the video with these parameters (Intensity `[0,200]`, Area `[150, 60000]`), it fails to pick up on one of the fishes during part of the video, and creates a "shadow" segment for one or both fishes.
This may be caused by the fact that neither fish moves for half of the video.
Try tracking with Intensity `[0,220]` which prevents the loss of segment for the second (MIKK) fish, but creates more "shadows".

View tracked video to confirm result.
