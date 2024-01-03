import os
import stat

from output_viewer.build import build_page, build_viewer
from output_viewer.index import (
    OutputFile,
    OutputGroup,
    OutputIndex,
    OutputPage,
    OutputRow,
)
from output_viewer.utils import rechmod


class OutputViewer(object):
    def __init__(self, path=".", index_name="Results"):
        self.path = os.path.abspath(path)
        self.index = OutputIndex(index_name)
        self.cache = {}  # dict of { OutputPage: { OutputGroup: [OutputRow] } }
        self.page = None
        self.group = None
        self.row = None

    def add_page(self, page_title, *args, **kwargs):
        """Add a page to the viewer's index"""
        self.page = OutputPage(page_title, *args, **kwargs)
        self.cache[self.page] = {}
        self.index.addPage(self.page)

    def set_page(self, page_title):
        """Sets the page with the title name as the current page"""
        for output_page in self.cache:
            if page_title == output_page.title:
                self.page = output_page
                return
        raise RuntimeError("There is no page titled: %s" % page_title)

    def add_group(self, group_name):
        """Add a group to the current page"""
        if self.page is None:
            raise RuntimeError("You must first insert a page with add_page()")
        self.group = OutputGroup(group_name)
        if self.group not in self.cache[self.page]:
            self.cache[self.page][self.group] = []  # group doesn't have any rows yet
        self.page.addGroup(self.group)

    def set_group(self, group_name):
        """Sets the group with the title name as the current group"""
        for output_group in self.cache[self.page]:
            if group_name == output_group.title:
                self.group = output_group
                return
        raise RuntimeError("There is no group titled: %s" % group_name)

    def add_row(self, row_name):
        """Add a row with the title name to the current group"""
        if self.group is None:
            raise RuntimeError("You must first insert a group with add_group()")
        self.row = OutputRow(row_name, [])
        if self.row not in self.cache[self.page][self.group]:
            self.cache[self.page][self.group].append(self.row)
        self.page.addRow(self.row, len(self.page.groups) - 1)  # type: ignore

    def set_row(self, row_name):
        """Sets the row with the title name as the current row"""
        for output_row in self.cache[self.page][self.group]:
            if row_name == output_row.title:
                self.row = output_row
                return
        raise RuntimeError("There is no row titled: %s" % row_name)

    def add_cols(self, cols):
        """Add multiple string cols to the current row"""
        self.row.columns.append(cols)  # type: ignore

    def add_col(self, col, is_file=False, **kwargs):
        """Add a single col to the current row. Set is_file to True if the col is a file path."""
        if is_file:
            self.row.columns.append(OutputFile(col, **kwargs))  # type: ignore
        else:
            self.row.columns.append(col)  # type: ignore

    def generate_page(self):
        """
        Generate and return the location of the current HTML page.
        """
        self.index.toJSON(os.path.join(self.path, "index.json"))

        default_mask = stat.S_IMODE(os.stat(self.path).st_mode)
        rechmod(self.path, default_mask)

        if os.access(self.path, os.W_OK):
            default_mask = stat.S_IMODE(
                os.stat(self.path).st_mode
            )  # mode of files to be included
            url = build_page(
                self.page,
                os.path.join(self.path, "index.json"),
                default_mask=default_mask,
            )
            return url

        raise RuntimeError("Error geneating the page.")

    def generate_viewer(self, prompt_user=True):
        """Generate the webpage and ask the user if they want to see it."""
        self.index.toJSON(os.path.join(self.path, "index.json"))

        default_mask = stat.S_IMODE(os.stat(self.path).st_mode)
        rechmod(self.path, default_mask)

        if os.access(self.path, os.W_OK):
            default_mask = stat.S_IMODE(
                os.stat(self.path).st_mode
            )  # mode of files to be included
            build_viewer(
                os.path.join(self.path, "index.json"),
                diag_name="My Diagnostics",
                default_mask=default_mask,
            )

        if os.path.exists(os.path.join(self.path, "index.html")):
            if prompt_user:
                user_prompt = "Would you like to open in a browser? y/[n]: "
                should_open = input(user_prompt)
                if should_open and should_open.lower()[0] == "y":
                    import webbrowser

                    index_path = os.path.join(self.path, "index.html")
                    url = "file://{path}".format(path=index_path)
                    webbrowser.open(url)
        else:
            raise RuntimeError("Failed to generate the viewer.")
