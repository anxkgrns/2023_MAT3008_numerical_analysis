import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

def expand_image(image, x_times ,y_times):
    height, width, channel = image.shape

    pixel = np.array(image)
    new_w = np.floor(x_times * width).astype(int)
    new_h = np.floor(y_times * height).astype(int)

    # 새로운 픽셀의 x,y
    new_pixel_x = np.array([[i for i in range(new_w)] for j in range(new_h)])
    new_pixel_y = np.array([[j for i in range(new_w)] for j in range(new_h)])


    # 원래 픽셀로 변환했을 때의 픽셀의 x,y
    original_pixel_x = (new_pixel_x/x_times)
    original_pixel_y = (new_pixel_y/y_times)
    
    # 원래 픽셀값의 가우스값
    gauss_x = np.floor(original_pixel_x).astype(int)
    gauss_y = np.floor(original_pixel_y).astype(int)
    
    gauss_x = np.minimum(gauss_x,width-2)
    gauss_y = np.minimum(gauss_y,height-2)
    

    # 가우스값과의 차이 - linear_interpolation.jpeg 참고
    diff_x = original_pixel_x - gauss_x
    diff_y = original_pixel_y - gauss_y

    # rgb 행렬로 스케일 맞추기 
    diff_x = np.repeat(np.expand_dims(diff_x, axis=-1), 3, axis=-1)
    diff_y = np.repeat(np.expand_dims(diff_y, axis=-1), 3, axis=-1)
    
    new_pixel = (1-diff_x)*(1-diff_y)*pixel[gauss_y,gauss_x] + (1-diff_x)*diff_y*pixel[gauss_y,gauss_x+1] + diff_x*(1-diff_y)*pixel[gauss_y+1,gauss_x] + diff_x*diff_y*pixel[gauss_y+1,gauss_x+1]
    return new_pixel


image = plt.imread("cat.jpg")
print("확대하고 싶은 배수를 입력하세요 ex) 2.1 1.5")
x_times, y_times = map(float,input().split()) 
new_image = expand_image(image, x_times, y_times)
Image.fromarray(new_image.astype(np.uint8)).save("expanded_image.jpg")

# 첫 번째 그래프 창
fig1, ax1 = plt.subplots(figsize=(image.shape[1]/100,image.shape[0]/100)) # (height, width) (이미지의 크기와 같은 크기로 그래프 창을 만들고 싶다면, figsize=
ax1.imshow(image)
ax1.set_title("Original Image")

# 두 번째 그래프 창
fig2, ax2 = plt.subplots(figsize=(new_image.shape[1]/100,new_image.shape[0]/100)) 
ax2.imshow(new_image.astype(np.uint8))

ax2.set_title("Expanded Image")

# 그래프 표시
plt.show()