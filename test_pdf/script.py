from pypdf import PdfWriter, PdfReader, PageObject

pdf_filenames = ["heatmap_alleles_alphas.pdf", "heatmap_populations_alphas.pdf"]

input1 = PdfReader(open(pdf_filenames[0], "rb"), strict=False)
input2 = PdfReader(open(pdf_filenames[1], "rb"), strict=False)

page1 = input1._get_page(0)
page2 = input2._get_page(0)

total_width = page1.mediabox.right + page2.mediabox.right
total_height = max([page1.mediabox.top, page2.mediabox.top])

new_page = PageObject.create_blank_page(None, total_width, total_height)

# Add first page at the 0,0 position
new_page.merge_page(page1)
# Add second page with moving along the axis x
new_page.merge_translated_page(page2, page1.mediabox.right, 0)

output = PdfWriter()
output.add_page(new_page)
output.write(open("result.pdf", "wb"))