
from django.template.base import Library

register = Library()

@register.filter("placeholder")
def add_placeholder(form, text=None):
    if text == None:
        raise ValueError("Placeholder requires text content for widget.")
    form.field.widget.attrs.update({ "placeholder": text })
    return form

@register.filter("autofocus")
def add_autofocus(form):
    form.field.widget.attrs.update({ "autofocus": "autofocus" })
    return form

@register.filter("tabindex")
def add_tabindex(form, number):
    form.field.widget.attrs.update({ "tabindex": number })
    return form

@register.filter("formfield")
def add_form_control(form):
    cls = ['form-control']
    if form.errors:
        cls.append("has-error")
    form.field.widget.attrs.update({"class": ' '.join(cls)})
    return form

