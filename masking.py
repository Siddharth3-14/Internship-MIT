from tkinter import *
from tkinter.filedialog import askopenfilename
from PIL import Image, ImageTk, ImageDraw

class Paint:
    def __init__(self, root):
        self.root = root
        self.root.title("Paint Program")
        self.brush_size = 3
        self.brush_color = "black"
        self.image_path = None

        # Create Canvas
        self.canvas = Canvas(root, bg="white", width=500, height=500)
        self.canvas.pack()

        # Create buttons
        self.clear_button = Button(root, text="Clear", command=self.clear_canvas)
        self.clear_button.pack(side=LEFT)

        self.choose_size_button = Scale(root, from_=1, to=10, orient=HORIZONTAL)
        self.choose_size_button.pack(side=LEFT)

        self.color_button = Button(root, text="Color", command=self.choose_color)
        self.color_button.pack(side=LEFT)

        self.open_image_button = Button(root, text="Open Image", command=self.open_image)
        self.open_image_button.pack(side=LEFT)

        # Bind mouse events to canvas
        self.canvas.bind("<B1-Motion>", self.draw)

    def draw(self, event):
        x1, y1 = (event.x - self.brush_size), (event.y - self.brush_size)
        x2, y2 = (event.x + self.brush_size), (event.y + self.brush_size)
        self.canvas.create_oval(x1, y1, x2, y2, fill=self.brush_color, outline=self.brush_color)

    def choose_color(self):
        self.brush_color = askcolor(color=self.brush_color)[1]

    def clear_canvas(self):
        self.canvas.delete("all")

    def open_image(self):
        self.image_path = askopenfilename(filetypes=[("Image Files", "*.png;*.jpg;*.jpeg")])
        image = Image.open(self.image_path)
        self.photo = ImageTk.PhotoImage(image)
        self.canvas.create_image(0, 0, image=self.photo, anchor=NW)

    def save_image(self):
        filename = asksaveasfilename(defaultextension=".png")
        self.canvas.postscript(file=filename, colormode="color")
        image = Image.open(filename)
        image.save(filename)

root = Tk()
paint_app = Paint(root)
root.mainloop()
