import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style('white')
# sns.set(font_scale=1.3)
plt.rcParams["font.family"] = "Times New Roman"
# Sample data
x = [1, 2, 3, 4, 5]
y = [2, 4, 6, 8, 10]

# Plotting the data
plt.plot(x, y, label='Line Plot')

# Adding a bold title 'A' outside the plot to the left
plt.text(0.5, 1.15, '(a)', transform=plt.gca().transAxes, fontsize=27, fontweight='bold', va='top', ha='right')

# Adding labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
# plt.title('Example Plot')

# Display the legend
plt.legend()
plt.savefig('plot.pdf', format='pdf', bbox_inches="tight")
# Show the plot
plt.show()