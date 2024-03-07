import PyPDF2


def add_header_footer_pdf(input_file, output_file, header_text):
    with open(input_file, 'rb') as file:
        pdf_reader = PyPDF2.PdfReader(file)
        pdf_writer = PyPDF2.PdfFileWriter()

        for page_num in range(pdf_reader.getNumPages()):
            page = pdf_reader.getPage(page_num)

            header = PyPDF2.pdf.PageObject.createBlankPage(None, page.mediaBox.getWidth(), 30)
            header.mergeTranslatedPage(page, 0, 30)
            header.mergeTranslatedPage(PyPDF2.pdf.PageObject.createTextObject(None, header_text), 0, 5)
            pdf_writer.addPage(header)

        with open(output_file, 'wb') as output:
            pdf_writer.write(output)


if __name__ == '__main__':
    add_header_footer_pdf('heatmap_alleles_alphas.pdf', 'Output.pdf', 'YasserKhalil')
